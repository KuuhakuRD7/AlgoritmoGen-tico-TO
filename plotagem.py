import matplotlib.pyplot as plt  # Importa a biblioteca para visualização de dados
import numpy as np  # Importa a biblioteca para manipulação de arrays e operações matemáticas

# Definindo a função F6, que calcula um valor baseado nas coordenadas x e y
def f6(x, y):
    return 0.5 - (np.sin(np.sqrt(x**2 + y**2))**2 - 0.5) / (1 + 0.001 * (x**2 + y**2)**2)

def plotar_F6_G(escolha):
    upperb = 10  # Limite superior para os eixos x e y
    lowerb = -10  # Limite inferior para os eixos x e y
    dados_por_geracao = {}  # Dicionário para armazenar dados por geração

    # Define o nome do arquivo com base na escolha do usuário
    if escolha == 0: 
        nome = "f6_geracoes_real.txt"  # Arquivo para dados de representação real
    elif escolha == 1: 
        nome = "f6_geracoes_bin.txt"  # Arquivo para dados de representação binária

    # Lê os dados do arquivo especificado
    with open(nome, "r") as file:
        for line in file:
            if line.strip():  # Verifica se a linha não está vazia
                geracao, x, y, fitness = map(float, line.split())  # Extrai os dados da linha
                if geracao not in dados_por_geracao:
                    dados_por_geracao[geracao] = ([], [], [])  # Inicializa listas para x, y e fitness
                dados_por_geracao[geracao][0].append(x)  # Adiciona x à lista correspondente
                dados_por_geracao[geracao][1].append(y)  # Adiciona y à lista correspondente
                dados_por_geracao[geracao][2].append(fitness)  # Adiciona fitness à lista correspondente

    # Criando uma grade para a função F6
    x_range = np.linspace(lowerb, upperb, 100)  # Gera 100 pontos entre lowerb e upperb para x
    y_range = np.linspace(lowerb, upperb, 100)  # Gera 100 pontos entre lowerb e upperb para y
    X, Y = np.meshgrid(x_range, y_range)  # Cria uma grade 2D a partir dos valores de x e y
    Z = f6(X, Y)  # Calcula os valores de Z usando a função F6

    # Plotando gráficos separados para cada geração
    for geracao, (x_values, y_values, fitness_values) in dados_por_geracao.items():
        fig = plt.figure(figsize=(12, 8))  # Cria uma nova figura com tamanho especificado
        ax = fig.add_subplot(111, projection='3d')  # Adiciona um subplot 3D

        # Plotando a superfície da função F6
        ax.plot_surface(X, Y, Z, cmap='coolwarm', alpha=0.5)  # Plota a superfície com um mapa de cores

        # Plotando os pontos da população para a geração atual
        ax.scatter(x_values, y_values, fitness_values, label=f'Geração {int(geracao)}', s=100, color='k', alpha=0.8)  # Plota os pontos da população

        # Configurando o gráfico
        ax.set_title(f'Gráfico 3D da Função F6 - Geração {int(geracao)}', fontsize=16)  # Título do gráfico
        ax.set_xlabel('X', fontsize=14)  # Rótulo do eixo X
        ax.set_ylabel('Y', fontsize=14)  # Rótulo do eixo Y
        ax.set_zlabel('Fitness', fontsize=14)  # Rótulo do eixo Z

        # Ajustando os limites dos eixos
        ax.set_xlim([lowerb, upperb])  # Limites do eixo X
        ax.set_ylim([lowerb, upperb])  # Limites do eixo Y
        ax.set_zlim([-1, 1])  # Limites do eixo Z (ajustar conforme necessário)

        # Ajustando a perspectiva do gráfico
        ax.view_init(elev=30, azim=30)  # Define a elevação e azimute para melhor visualização

        ax.legend()  # Adiciona a legenda para identificar a geração
        # # Salva o gráfico como um arquivo PNG (opcional)
        # if escolha == 1:
        #     plt.savefig(f"bin_{int(geracao)}.png", bbox_inches='tight')  # Salva a figura para binário
        # elif escolha == 0:
        #     plt.savefig(f"real_{int(geracao)}.png", bbox_inches='tight')  # Salva a figura para real

        plt.show()  # Exibe o gráfico

def plotar_indiv_G(data):
    # Extrai os valores de geração, F6, x e y
    generations = data[:, 0].astype(int)  # Primeira coluna: gerações
    F6_values = data[:, 1]  # Segunda coluna: valores de F6
    x_values = data[:, 2]    # Terceira coluna: valores de x
    y_values = data[:, 3]    # Quarta coluna: valores de y

    plt.figure(figsize=(20, 10))  # Aumenta o tamanho da figura
    best_fitness_per_generation = []  # Lista para armazenar o melhor fitness por geração

    for generation in range(len(F6_values)):  # Itera sobre as gerações
        best_fitness = F6_values[generation]  # Encontra o melhor fitness na geração atual
        best_fitness_per_generation.append(best_fitness)  # Adiciona o melhor fitness à lista

    # Plota a curva do melhor fitness
    plt.plot(best_fitness_per_generation, color='b', marker='o', markersize=3)  # Plota a curva do melhor fitness
    plt.title('Curva do Melhor Indivíduo Geração a Geração', fontsize=18)  # Título da plotagem
    plt.xlabel('Gerações', fontsize=16)  # Rótulo do eixo X
    plt.ylabel('Melhor Fitness', fontsize=16)  # Rótulo do eixo Y
    plt.grid()  # Adiciona uma grade
    plt.xlim(0, len(best_fitness_per_generation))  # Ajusta os limites do eixo X
    plt.ylim(min(best_fitness_per_generation) - 0.1, max(best_fitness_per_generation) + 0.1)  # Ajusta os limites do eixo Y

    plt.show()  # Exibe a plotagem da curva

def plotar_indiv_partes(data):
    # Extrai os valores de geração, F6, x e y
    generations = data[:, 0].astype(int)  # Primeira coluna: gerações
    F6_values = data[:, 1]  # Segunda coluna: valores de F6
    x_values = data[:, 2]    # Terceira coluna: valores de x
    y_values = data[:, 3]    # Quarta coluna: valores de y

    # Define os intervalos de gerações para plotagem
    intervals = [(0, 90), (100, 190), (200, 290), (300, 390), (400, 490)]
    titles = ['Gerações 0-90', 'Gerações 100-190', 'Gerações 200-290', 'Gerações 300-390', 'Gerações 400-490']

    for i, (start, end) in enumerate(intervals):
        plt.figure(figsize=(20, 10))  # Cria uma nova figura para cada intervalo
        best_fitness_per_generation = []  # Lista para armazenar o melhor fitness por geração

        for generation in range(start, end + 1):
            if generation < len(F6_values):  # Verifica se a geração está dentro do limite
                best_fitness = F6_values[generation]  # Encontra o melhor fitness na geração atual
                best_fitness_per_generation.append(best_fitness)  # Adiciona o melhor fitness à lista

        # Plota a curva do melhor fitness
        plt.plot(best_fitness_per_generation, color='b', marker='o', markersize=3)  # Plota a curva do melhor fitness
        plt.title(titles[i], fontsize=18)  # Título da plotagem
        plt.xlabel('Gerações', fontsize=16)  # Rótulo do eixo X
        plt.ylabel('Melhor Fitness', fontsize=16)  # Rótulo do eixo Y
        plt.grid()  # Adiciona uma grade
        plt.xlim(0, len(best_fitness_per_generation))  # Ajusta os limites do eixo X
        plt.ylim(min(best_fitness_per_generation) - 0.1, max(best_fitness_per_generation) + 0.1)  # Ajusta os limites do eixo Y

        plt.show()  # Exibe a plotagem da curva para o intervalo atual

def plotar_matrix():
    # Lê os dados do arquivo que contém as matrizes
    data = np.genfromtxt("matrizes.txt", dtype=str, delimiter=' ', invalid_raise=False)
    # Remove linhas vazias e mantém apenas as que contêm valores binários
    data = data[~np.char.equal(data[:, 0], '')]  # Remove linhas vazias

    # Converte a primeira coluna para inteiros (gerações)
    data[:, 0] = data[:, 0].astype(int)

    # Obtém as representações binárias de X e Y
    geracoes = np.unique(data[:, 0])  # Obtemos as gerações únicas
    for geracao in geracoes:
        # Filtra os dados para a geração atual
        dados_geracao = data[data[:, 0] == geracao]
        x_bin = dados_geracao[:, 1]  # Primeira coluna: valores binários de X

        # Converte para matriz binária para exibição
        binary_matrix = np.array([[int(bit) for bit in list(row)] for row in x_bin[:200]])  # 200 indivíduos

        # Cria um gráfico usando imshow
        plt.figure(figsize=(20, 10))  # Aumenta o tamanho da figura
        plt.imshow(binary_matrix, cmap='gray', aspect='auto')  # Usa cmap='gray' para representar 0s e 1s
        plt.title(f"Gráfico de Dispersão da População - Geração {int(geracao)}", fontsize=18)
        plt.xlabel("Bits do Cromossomo", fontsize=16)  # Rótulo do eixo X
        plt.ylabel("Indivíduos", fontsize=16)  # Rótulo do eixo Y
        plt.colorbar(label='Valor Binário')  # Adiciona uma barra de cores
        plt.xlim(-1, binary_matrix.shape[1])  # Ajusta os limites do eixo X
        plt.ylim(-1, binary_matrix.shape[0])  # Ajusta os limites do eixo Y
        plt.gca().invert_yaxis()  # Inverte o eixo Y para que o primeiro indivíduo fique no topo
        plt.grid(False)  # Desativa a grade para melhor visualização

        # Salva o gráfico como um arquivo PNG
        plt.savefig(f"matriz_{int(geracao)}.png", bbox_inches='tight')  # Salva a figura
        plt.close()  # Fecha a figura para liberar memória

# Chamada das funções principais
escolha = int(input("Escolha 0 para plotar o AG - Real e 1 para plotar o AG - Binário: \n"))
if escolha == 0:
    data = np.loadtxt('results_real.txt')  # Carrega os dados para representação real
elif escolha == 1:
    data = np.loadtxt('results_bin.txt')  # Carrega os dados para representação binária
else:
    print("Escolha Inválida\n")  # Mensagem de erro para escolha inválida
    exit()


# Chamada das funções de plotagem (descomente conforme necessário uma por vez)
# plotar_indiv_G(data)
# plotar_indiv_partes(data)
# plotar_F6_G(escolha)
# if escolha == 1: plotar_matrix()
