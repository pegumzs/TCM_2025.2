import numpy as np
import matplotlib.pyplot as plt

# Configurações globais de plotagem
plt.rcParams['figure.figsize'] = [10, 6]
plt.rcParams['font.size'] = 12
plt.rcParams['axes.grid'] = True

def solucao_analitica(x, t, problema, n_termos=100):
    """Calcula a solução exata ."""
    if problema == 1:

        soma = np.zeros_like(x)
        for n in range(1, n_termos + 1):
            termo = (4 / ((2*n - 1) * np.pi)) * \
                    np.sin((2*n - 1) * np.pi * x) * \
                    np.exp(-((2*n - 1)**2) * (np.pi**2) * t)
            soma += termo
        return soma

    elif problema == 2:

        soma = np.zeros_like(x)
        for n in range(1, n_termos + 1):
            termo = (2 / (n * np.pi)) * \
                    np.sin(n * np.pi * x) * \
                    np.exp(-(n**2) * (np.pi**2) * t)
            soma += termo
        return 1 - x - soma

    elif problema == 3:

        return np.exp(-(np.pi**2 * t) / 4) * np.sin((np.pi * x) / 2)

def resolver_problema(id_problema, L, N, alpha, tempos_analise, titulo):
    """
    Resolve a equação do calor 1D numérico vs analítico.
    """

    # 1. Configuração da Malha e Parâmetros
    dx = L / N
    # Critério de estabilidada
    dt = 0.2 * (dx**2) / alpha
    lamb = alpha * dt / (dx**2)  # Número de Fourier

    print(f"L={L}, N={N}, dx={dx:.4f}")
    print(f"dt={dt:.5f}, Lambda={lamb:.4f} (Estável se <= 0.5)")

    # Malha espacial
    x = np.linspace(0, L, N + 1)

    # 2. Condições Iniciais (t=0)
    T = np.zeros(N + 1)

    if id_problema == 1:
        T[:] = 1.0        # Barra quente
        T[0], T[-1] = 0.0, 0.0 # Bordas frias
    elif id_problema == 2:
        T[:] = 0.0        # Barra fria
        T[0] = 1.0        # Borda esq quente
        T[-1] = 0.0       # Borda dir fria
    elif id_problema == 3:
        T = np.sin((np.pi * x) / 2) # Senoide
        T[0], T[-1] = 0.0, 0.0

    # Cópias para o loop
    T_atual = np.copy(T)
    T_nova = np.copy(T)

    # Prepara o gráfico
    plt.figure()
    plt.title(titulo)
    plt.xlabel("Posição x")
    plt.ylabel("Temperatura T")

    # Cores para diferenciar os tempos
    cores = ['b', 'g', 'r', 'm', 'k', 'c']

    # Loop principal de tempo
    tempo_simulado = 0.0
    max_tempo = max(tempos_analise)

    # Índice para controlar as cores e labels
    idx_plot = 0

    # Usando WHILE para controle contínuo do tempo (recomendação Video )
    while tempo_simulado < max_tempo:

        for i in range(1, N):
            T_nova[i] = T_atual[i] + lamb * (T_atual[i+1] - 2*T_atual[i] + T_atual[i-1])

        # Atualiza vetor e tempo
        T_atual = np.copy(T_nova)
        tempo_simulado += dt

        # Verifica se atingiu um tempo de análise
        # (Usamos uma lista ordenada para plotar na ordem)
        if idx_plot < len(tempos_analise):
            t_alvo = tempos_analise[idx_plot]


            if tempo_simulado >= t_alvo:
                cor = cores[idx_plot % len(cores)]

                # Plot Numérico
                plt.plot(x, T_atual, 'o', color=cor, markersize=5,
                         label=f'Num t={t_alvo:.2f}')

                # Plot Analítico
                T_exata = solucao_analitica(x, t_alvo, id_problema)
                plt.plot(x, T_exata, '-', color=cor, linewidth=1.5, alpha=0.7)

                idx_plot += 1

    plt.plot([], [], 'ko', label='Numérico')
    plt.plot([], [], 'k-', label='Analítico')

    plt.legend(loc='best', fontsize='small', ncol=2)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.show()


def resolver_numerico(id_problema, L, N, alpha, tempo_alvo):
    # 1. Configuração
    dx = L / N
    # Mantendo o mesmo critério de estabilidade do código original
    dt = 0.2 * (dx**2) / alpha
    lamb = alpha * dt / (dx**2)

    # 2. Malha e Condições Iniciais
    x = np.linspace(0, L, N + 1)
    T = np.zeros(N + 1)

    # Configuração para o Problema 1 (usado na análise de Alpha)
    if id_problema == 1:
        T[:] = 1.0
        T[0], T[-1] = 0.0, 0.0
    # (Adicione elif para outros problemas se for variar o problema na análise)

    # 3. Loop no Tempo
    T_atual = np.copy(T)
    T_nova = np.copy(T)
    tempo_simulado = 0.0

    while tempo_simulado < tempo_alvo:
        for i in range(1, N):
            T_nova[i] = T_atual[i] + lamb * (T_atual[i+1] - 2*T_atual[i] + T_atual[i-1])

        T_atual = np.copy(T_nova)
        tempo_simulado += dt

    return x, T_atual

# PROBLEMA 1: Resfriamento
# Tempos estendidos para garantir visualização do regime permanente (t=0.5 -> T=0)
resolver_problema(
    id_problema=1,
    L=1.0, N=20, alpha=1.0,
    tempos_analise=[0.01, 0.05, 0.1, 0.5],
    titulo="Problema 1: Resfriamento (Bordas 0, Inicial 1)"
)

# PROBLEMA 2: Aquecimento Assimétrico
# Regime permanente é uma reta decrescente de 1 a 0
resolver_problema(
    id_problema=2,
    L=1.0, N=20, alpha=1.0,
    tempos_analise=[0.02, 0.1, 0.2, 0.5],
    titulo="Problema 2: Aquecimento Esq (Esq 1, Dir 0)"
)

# PROBLEMA 3: Senoide (L=2)
# Regime permanente é zero (decaimento exponencial)
resolver_problema(
    id_problema=3,
    L=2.0, N=20, alpha=1.0,
    tempos_analise=[0.1, 0.5, 1.0, 1.5],
    titulo="Problema 3: Senoide (L=2, Bordas 0)"
)


# 4.2 Influência de Alpha (Difusividade Térmica)

plt.figure()
t_fixo = 0.1
N_fixo = 51
alphas = [0.1, 1.0, 5.0] # Alpha baixo (lento), médio e alto (rápido)

for a_val in alphas:
    x_res, T_res = resolver_numerico(1, 1.0, N_fixo, a_val, t_fixo)
    plt.plot(x_res, T_res, linewidth=2, label=f'Alpha = {a_val}')

plt.title(f"Influência da Difusividade Térmica (t={t_fixo}s)")
plt.xlabel("Posição x")
plt.ylabel("Temperatura T")
plt.legend()
plt.savefig("analise_alpha.png", dpi=300)
print("Salvo: analise_alpha.png")
plt.show()
