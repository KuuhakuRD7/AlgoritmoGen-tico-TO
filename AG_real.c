#include <stdio.h>   // Biblioteca para entrada e saída padrão
#include <stdlib.h>  // Biblioteca para funções de alocação de memória e geração de números aleatórios
#include <math.h>    // Biblioteca para funções matemáticas
#include <time.h>    // Biblioteca para manipulação de tempo


// Estrutura que representa uma solução no algoritmo genético
typedef struct {
    double x; // Valor real de x
    double y; // Valor real de y
    double fitness_value; // Valor da função F6
} Solution;

// Função para calcular a função de fitness
void calculate_fitness(Solution* s) {
    // Calcula o valor da função F6
    s->fitness_value = 0.5 -( (pow(sin(sqrt(s->x * s->x + s->y * s->y)), 2) - 0.5) /
                          pow(1 + 0.001 * (s->x * s->x + s->y * s->y), 2) );
}

// Função para seleção por roleta
Solution roulette_selection(Solution* population, int population_size) {
    double total_fitness = 0.0; // Inicializa a soma total da aptidão
    // Calcula a soma total da aptidão da população
    for (int i = 0; i < population_size; i++) {
        total_fitness += population[i].fitness_value;
    }

    // Gera um valor aleatório proporcional à aptidão total
    double random_value = ((double)rand() / RAND_MAX) * total_fitness;
    double cumulative_sum = 0.0; // Inicializa a soma cumulativa

    // Seleciona um indivíduo com base na roleta
    for (int i = 0; i < population_size; i++) {
        cumulative_sum += population[i].fitness_value; // Atualiza a soma cumulativa
        if (cumulative_sum >= random_value) { // Se a soma cumulativa ultrapassar o valor aleatório
            return population[i]; // Retorna o indivíduo selecionado
        }
    }
    return population[population_size - 1]; // Retorna o último caso não encontrado
}

// // Função para cruzamento de um ponto
// void crossover_one_point(Solution parent1, Solution parent2, Solution* child1, Solution* child2, double Tc) {
//     if (((double)rand() / RAND_MAX) < Tc) { // Verifica se deve ocorrer o cruzamento com base na taxa de cruzamento
//         // Realiza o cruzamento para x e y
//         child1->x = (parent1.x + parent2.x) / 2; // Média dos valores de x
//         child1->y = (parent1.y + parent2.y) / 2; // Média dos valores de y
//
//         child2->x = (parent1.x + parent2.x) / 2; // Média dos valores de x
//         child2->y = (parent1.y + parent2.y) / 2; // Média dos valores de y
//     } else {
//         *child1 = parent1; // Se não houver cruzamento, o filho é igual ao pai
//         *child2 = parent2; // O mesmo para o segundo filho
//     }
// }
//
// // Função para cruzamento de dois pontos
// void crossover_two_points(Solution parent1, Solution parent2, Solution* child1, Solution* child2, double Tc) {
//     if (((double)rand() / RAND_MAX) < Tc) { // Verifica se deve ocorrer o cruzamento com base na taxa de cruzamento
//         // Escolhe dois pontos de corte aleatórios
//         double alpha = ((double)rand() / RAND_MAX);
//         double beta = ((double)rand() / RAND_MAX);
//
//         // Realiza o cruzamento para x e y
//         child1->x = alpha * parent1.x + (1 - alpha) * parent2.x; // Combina os valores de x
//         child1->y = alpha * parent1.y + (1 - alpha) * parent2.y; // Combina os valores de y
//
//         child2->x = beta * parent1.x + (1 - beta) * parent2.x; // Combina os valores de x
//         child2->y = beta * parent1.y + (1 - beta) * parent2.y; // Combina os valores de y
//     } else {
//         *child1 = parent1; // Se não houver cruzamento, o filho é igual ao pai
//         *child2 = parent2; // O mesmo para o segundo filho
//     }
// }
//
// // Função para cruzamento uniforme
// void crossover_uniform(Solution parent1, Solution parent2, Solution* child1, Solution* child2, double Tc) {
//     if (((double)rand() / RAND_MAX) < Tc) { // Verifica se deve ocorrer o cruzamento com base na taxa de cruzamento
//         // Realiza o cruzamento uniforme
//         child1->x = ((double)rand() / RAND_MAX < 0.5) ? parent1.x : parent2.x; // Escolhe aleatoriamente entre os pais
//         child1->y = ((double)rand() / RAND_MAX < 0.5) ? parent1.y : parent2.y; // Escolhe aleatoriamente entre os pais
//
//         child2->x = ((double)rand() / RAND_MAX < 0.5) ? parent1.x : parent2.x; // Escolhe aleatoriamente entre os pais
//         child2->y = ((double)rand() / RAND_MAX < 0.5) ? parent1.y : parent2.y; // Escolhe aleatoriamente entre os pais
//     } else {
//         *child1 = parent1; // Se não houver cruzamento, o filho é igual ao pai
//         *child2 = parent2; // O mesmo para o segundo filho
//     }
// }
//
// // Função para cruzamento aritmético
// void crossover_arithmetic_weighted(Solution parent1, Solution parent2, Solution* child1, Solution* child2, double Tc) {
//     if (((double)rand() / RAND_MAX) < Tc) { // Verifica se deve ocorrer o cruzamento com base na taxa de cruzamento
//         child1->x = 0.5 * (parent1.x + parent2.x); // Média aritmética dos valores de x
//         child1->y = 0.5 * (parent1.y + parent2.y); // Média aritmética dos valores de y
//
//         child2->x = 0.5 * (parent1.x + parent2.x); // Média aritmética dos valores de x
//         child2->y = 0.5 * (parent1.y + parent2.y); // Média aritmética dos valores de y
//     } else {
//         *child1 = parent1; // Se não houver cruzamento, o filho é igual ao pai
//         *child2 = parent2; // O mesmo para o segundo filho
//     }
// }

// Função para cruzamento por média aritmética
void crossover_arithmetic(Solution parent1, Solution parent2, Solution* child1, Solution* child2, double Tc) {
    if (((double)rand() / RAND_MAX) < Tc) { // Verifica se deve ocorrer o cruzamento com base na taxa de cruzamento
        // double alpha = ((double)rand() / RAND_MAX); // Peso aleatório
        double alpha = 0.35;
        child1->x = alpha * parent1.x + (1 - alpha) * parent2.x; // Combina os valores de x
        child1->y = alpha * parent1.y + (1 - alpha) * parent2.y; // Combina os valores de y

        child2->x = (1 - alpha) * parent1.x + alpha * parent2.x; // Combina os valores de x
        child2->y = (1 - alpha) * parent1.y + alpha * parent2.y; // Combina os valores de y
    } else {
        *child1 = parent1; // Se não houver cruzamento, o filho é igual ao pai
        *child2 = parent2; // O mesmo para o segundo filho
    }
}

// // Função para mutação randômica uniforme
// void mutate_uniform(Solution* s, double Tm, int generation, int max_generations) {
//     // Verifica se deve ocorrer mutação com base na taxa de mutação
//     if (((double)rand() / RAND_MAX) < Tm) { // Taxa de mutação de 1%
//         s->x += ((double)rand() / RAND_MAX - 0.5) * 2; // Aplica uma pequena perturbação em x
//         if (s->x < -100.0) s->x = -100.0;
//         if (s->x > 100.0) s->x = 100.0;
//     }
//     if (((double)rand() / RAND_MAX) < Tm) { // Taxa de mutação de 1%
//         s->y += ((double)rand() / RAND_MAX - 0.5) * 2; // Aplica uma pequena perturbação em y
//         if (s->y < -100.0) s->y = -100.0;
//         if (s->y > 100.0) s->y = 100.0;
//     }
// }

double gaussiano(double media, double variancia) {
    double u1 = (double)rand() / RAND_MAX;
    double u2 = (double)rand() / RAND_MAX;
    return media + sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2) * sqrt(variancia);
}

// Função para mutação randômica não uniforme
void mutate_non_uniform(Solution* s, double Tm, int generation, int max_generations) {
    double varinicial = 10;
    double variancia = varinicial * pow(1 - (double)generation / max_generations, 2);
    // Verifica se deve ocorrer mutação com base na taxa de mutação
    if (((double)rand() / RAND_MAX) < Tm * (1 - (double)generation / max_generations)) { // Taxa de mutação decrescente
        //s->x += ((double)rand() / RAND_MAX - 0.5) * 2; // Aplica uma pequena perturbação em x
        s->x += gaussiano(0, variancia);
        if (s->x < -100.0) s->x = -100.0;
        if (s->x > 100.0) s->x = 100.0;
    }
    if (((double)rand() / RAND_MAX) < Tm * (1 - (double)generation / max_generations)) { // Taxa de mutação decrescente
        //s->y += ((double)rand() / RAND_MAX - 0.5) * 2; // Aplica uma pequena perturbação em y
        s->y += gaussiano(0, variancia);
        if (s->y < -100.0) s->y = -100.0;
        if (s->y > 100.0) s->y = 100.0;
    }
}


// Função que executa o algoritmo genético
void run_genetic_algorithm(int population_size, int generations, const char* filename, double Tc, double Tm, void (*crossover_func)(Solution, Solution, Solution*, Solution*, double), void (*mutate_func)(Solution*, double, int, int)) {
    // Aloca memória para a população
    Solution* population = (Solution*)malloc(population_size * sizeof(Solution));
    if (population == NULL) { // Verifica se a alocação foi bem-sucedida
        fprintf(stderr, "Erro ao alocar memória.\n");
        return; // Sai da função em caso de erro
    }

    srand((unsigned int)time(NULL)); // Inicializa o gerador de números aleatórios

    // Inicializa a população com valores reais aleatórios
    for (int i = 0; i < population_size; i++) {
        population[i].x = ((double)rand() / RAND_MAX) * 200.0 - 100.0; // Gera um valor real aleatório para x
        population[i].y = ((double)rand() / RAND_MAX) * 200.0 - 100.0; // Gera um valor real aleatório para y
        calculate_fitness(&population[i]); // Calcula a aptidão da solução
    }

    // Cria um arquivo para salvar os resultados
    FILE* output_file = fopen(filename, "w");
    if (!output_file) { // Verifica se o arquivo foi aberto corretamente
        fprintf(stderr, "Erro ao abrir o arquivo para escrita.\n");
        free(population); // Libera a memória da população
        return; // Sai da função em caso de erro
    }

    // Variáveis para rastrear o melhor fitness
    double best_fitness = -INFINITY; // Inicializa com o menor valor possível
    double best_x = 0.0;
    double best_y = 0.0;
    int best_generation = 0;

    // Salva os valores de todos os indivíduos a cada 10 gerações
    FILE* f6_file = fopen("f6_geracoes_real.txt", "a");
    if (!f6_file) {
        fprintf(stderr, "Erro ao abrir o arquivo f6_geracoes.txt para escrita.\n");
        fclose(output_file);
        free(population);
        return;
    }

    // Loop para cada geração
    for (int generation = 0; generation < generations; ++generation) {
        // Aloca memória para a nova população
        Solution* new_population = (Solution*)malloc(population_size * sizeof(Solution));
        if (new_population == NULL) { // Verifica se a alocação foi bem-sucedida
            fprintf(stderr, "Erro ao alocar memória para nova população.\n");
            fclose(output_file); // Fecha o arquivo
            fclose(f6_file);
            free(population); // Libera a memória da população
            return; // Sai da função em caso de erro
        }

        // Gera novos indivíduos
        for (int i = 0; i < population_size; i += 2) {
            Solution parent1 = roulette_selection(population, population_size); // Seleciona o primeiro pai
            Solution parent2 = roulette_selection(population, population_size); // Seleciona o segundo pai
            Solution child1, child2; // Declara os filhos

            crossover_func(parent1, parent2, &child1, &child2, Tc); // Realiza o cruzamento usando a função escolhida
            mutate_func(&child1, Tm, generation, generations); // Aplica mutação ao primeiro filho
            mutate_func(&child2, Tm, generation, generations); // Aplica mutação ao segundo filho

            calculate_fitness(&child1); // Calcula a aptidão do primeiro filho
            calculate_fitness(&child2); // Calcula a aptidão do segundo filho

            new_population[i] = child1; // Adiciona o primeiro filho à nova população
            new_population[i + 1] = child2; // Adiciona o segundo filho à nova população
        }

        // Atualiza a população com a nova população
        for (int i = 0; i < population_size; i++) {
            population[i] = new_population[i];
        }


        // Salva os valores de todos os indivíduos a cada 10 gerações
        if (generation % 10 == 0) {
            for (int j = 0; j < population_size; j++) {
                fprintf(f6_file, "%d %f %f %f\n", generation, population[j].x, population[j].y, population[j].fitness_value);
            }
        }

        // Encontra a melhor solução da geração atual
        Solution best_solution = population[0]; // Assume que o primeiro é o melhor
        for (int i = 1; i < population_size; i++) {
            if (population[i].fitness_value > best_solution.fitness_value) { // Se encontrar uma melhor
                best_solution = population[i]; // Atualiza a melhor solução
            }
        }

        // Atualiza o melhor fitness encontrado
        if (best_solution.fitness_value > best_fitness) {
            best_fitness = best_solution.fitness_value;
            best_x = best_solution.x;
            best_y = best_solution.y;
            best_generation = generation; // Armazena a geração do melhor fitness
        }

        // Salva os resultados no arquivo
        fprintf(output_file, "%d %f %f %f\n", generation, best_solution.fitness_value, best_solution.x, best_solution.y);
        // Imprime a melhor solução da geração atual
        printf("Geração %d: Melhor Fitness = %f | x = %f | y = %f\n", generation, best_solution.fitness_value, best_solution.x, best_solution.y);

        free(new_population); // Libera a memória da nova população
    }

    // // Salva o melhor fitness encontrado e suas coordenadas
    // fprintf(output_file, "%d %f %f %f\n",
    //         best_generation, best_fitness, best_x, best_y);
    //
    printf("Melhor Fitness encontrado: %f em Geração %d | x = %f | y = %f\n",
             best_fitness, best_generation, best_x, best_y);

    fclose(output_file); // Fecha o arquivo
    fclose(f6_file); // Fecha o arquivo após a escrita
    free(population); // Libera a memória da população
}

int main() {
    // Chama a função do algoritmo genético com parâmetros específicos e o método de cruzamento e mutação desejados
    //run_genetic_algorithm(num_individuos, num_geracoes, taxa_cruzamento, taxa_mutacao, tipo_cruzamento, tipo_mutacao)
    run_genetic_algorithm(200, 500, "results_real.txt", 0.85, 0.01, crossover_arithmetic, mutate_non_uniform); // Executa com cruzamento aritmético e mutação uniforme

    return 0; // Retorna 0 para indicar que o programa terminou com sucesso
}
