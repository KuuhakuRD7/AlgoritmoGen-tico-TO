#include <stdio.h>   // Biblioteca para entrada e saída padrão
#include <stdlib.h>  // Biblioteca para funções de alocação de memória e geração de números aleatórios
#include <math.h>    // Biblioteca para funções matemáticas
#include <time.h>    // Biblioteca para manipulação de tempo

#define BITS 25 // Define o número de bits para a representação binária

// Estrutura que representa uma solução no algoritmo genético
typedef struct {
    unsigned int x_bin; // Representação binária de x
    unsigned int y_bin; // Representação binária de y
    double fitness_value; // Valor da função F6
} Solution;

// Função para decodificar o valor binário para o intervalo [-100, 100]
double decode(unsigned int value) {
    return (value / (double)((1 << BITS) - 1)) * 200.0 - 100.0; // Mapeia o valor binário para o intervalo [-100, 100]
}

// Função para converter um número inteiro em uma string de sua representação binária
void int_to_binary_string(unsigned int value, char* buffer, int buffer_size) {
    for (int i = 0; i < buffer_size - 1; i++) {
        buffer[buffer_size - 2 - i] = (value & (1 << i)) ? '1' : '0';
    }
    buffer[buffer_size - 1] = '\0'; // Null-terminate the string
}

// Função para calcular a função de fitness
void calculate_fitness(Solution* s) {
    double x = decode(s->x_bin); // Decodifica x
    double y = decode(s->y_bin); // Decodifica y
    // Calcula o valor da função F6
    s->fitness_value = 0.5 - (pow(sin(sqrt(x * x + y * y)), 2) - 0.5) /
                          pow(1 + 0.001 * (x * x + y * y), 2);
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

// Função para cruzamento de um ponto
void crossover_one_point(Solution parent1, Solution parent2, Solution* child1, Solution* child2, double Tc) {
    if (((double)rand() / RAND_MAX) < Tc) { // Verifica se deve ocorrer o cruzamento com base na taxa de cruzamento
        int crossover_point = rand() % BITS; // Ponto de cruzamento em 25 bits

        // Realiza o cruzamento para x
        child1->x_bin = (parent1.x_bin & ((1 << crossover_point) - 1)) | (parent2.x_bin & ~((1 << crossover_point) - 1));
        // Realiza o cruzamento para y
        child1->y_bin = (parent1.y_bin & ((1 << crossover_point) - 1)) | (parent2.y_bin & ~((1 << crossover_point) - 1));

        // Realiza o cruzamento para o segundo filho
        child2->x_bin = (parent2.x_bin & ((1 << crossover_point) - 1)) | (parent1.x_bin & ~((1 << crossover_point) - 1));
        child2->y_bin = (parent2.y_bin & ((1 << crossover_point) - 1)) | (parent1.y_bin & ~((1 << crossover_point) - 1));
    } else {
        *child1 = parent1; // Se não houver cruzamento, o filho é igual ao pai
        *child2 = parent2; // O mesmo para o segundo filho
    }
}

// // Função para cruzamento de dois pontos
// void crossover_two_points(Solution parent1, Solution parent2, Solution* child1, Solution* child2, double Tc) {
//     if (((double)rand() / RAND_MAX) < Tc) { // Verifica se deve ocorrer o cruzamento com base na taxa de cruzamento
//         int point1 = rand() % BITS; // Primeiro ponto de cruzamento
//         int point2 = rand() % BITS; // Segundo ponto de cruzamento
//         if (point1 > point2) { // Garante que point1 seja menor que point2
//             int temp = point1;
//             point1 = point2;
//             point2 = temp;
//         }
//
//         // Realiza o cruzamento para x
//         child1->x_bin = (parent1.x_bin & ((1 << point1) - 1)) | (parent2.x_bin & ~((1 << point1) - 1));
//         child1->x_bin |= (parent2.x_bin & ((1 << (point2 + 1)) - 1)) | (parent1.x_bin & ~((1 << (point2 + 1)) - 1));
//
//         // Realiza o cruzamento para y
//         child1->y_bin = (parent1.y_bin & ((1 << point1) - 1)) | (parent2.y_bin & ~((1 << point1) - 1));
//         child1->y_bin |= (parent2.y_bin & ((1 << (point2 + 1)) - 1)) | (parent1.y_bin & ~((1 << (point2 + 1)) - 1));
//
//         // Realiza o cruzamento para o segundo filho
//         child2->x_bin = (parent2.x_bin & ((1 << point1) - 1)) | (parent1.x_bin & ~((1 << point1) - 1));
//         child2->x_bin |= (parent1.x_bin & ((1 << (point2 + 1)) - 1)) | (parent2.x_bin & ~((1 << (point2 + 1)) - 1));
//
//         // Realiza o cruzamento para y
//         child2->y_bin = (parent2.y_bin & ((1 << point1) - 1)) | (parent1.y_bin & ~((1 << point1) - 1));
//         child2->y_bin |= (parent1.y_bin & ((1 << (point2 + 1)) - 1)) | (parent2.y_bin & ~((1 << (point2 + 1)) - 1));
//     } else {
//         *child1 = parent1; // Se não houver cruzamento, o filho é igual ao pai
//         *child2 = parent2; // O mesmo para o segundo filho
//     }
// }
//
// // Função para cruzamento uniforme
// void crossover_uniform(Solution parent1, Solution parent2, Solution* child1, Solution* child2, double Tc) {
//     if (((double)rand() / RAND_MAX) < Tc) { // Verifica se deve ocorrer o cruzamento com base na taxa de cruzamento
//         child1->x_bin = 0; // Inicializa child1
//         child1->y_bin = 0; // Inicializa child1
//         child2->x_bin = 0; // Inicializa child2
//         child2->y_bin = 0; // Inicializa child2
//
//         for (int i = 0; i < BITS; i++) {
//             if (rand() % 2) {
//                 child1->x_bin |= (parent1.x_bin & (1 << i)); // Escolhe aleatoriamente entre os pais
//                 child1->y_bin |= (parent1.y_bin & (1 << i)); // Escolhe aleatoriamente entre os pais
//             } else {
//                 child1->x_bin |= (parent2.x_bin & (1 << i));
//                 child1->y_bin |= (parent2.y_bin & (1 << i));
//             }
//         }
//         // O mesmo para o segundo filho
//         for (int i = 0; i < BITS; i++) {
//             if (rand() % 2) {
//                 child2->x_bin |= (parent1.x_bin & (1 << i));
//                 child2->y_bin |= (parent1.y_bin & (1 << i));
//             } else {
//                 child2->x_bin |= (parent2.x_bin & (1 << i));
//                 child2->y_bin |= (parent2.y_bin & (1 << i));
//             }
//         }
//     } else {
//         *child1 = parent1; // Se não houver cruzamento, o filho é igual ao pai
//         *child2 = parent2; // O mesmo para o segundo filho
//     }
// }

// Função para mutação
void mutate(Solution* s, double Tm) {
    for (int i = 0; i < BITS; ++i) {
        // Verifica se deve ocorrer mutação com base na taxa de mutação
        if (((double)rand() / RAND_MAX) < Tm) { // Taxa de mutação
            s->x_bin ^= (1 << i); // Inverte o bit de x
        }
        if (((double)rand() / RAND_MAX) < Tm) { // Taxa de mutação
            s->y_bin ^= (1 << i); // Inverte o bit de y
        }
    }
}

// // Função para mutação randômica uniforme
// void mutate_uniform(Solution* s, double Tm) {
//     if (((double)rand() / RAND_MAX) < Tm) {
//         s->x_bin ^= (1 << (rand() % BITS)); // Inverte um bit aleatório de x
//     }
//     if (((double)rand() / RAND_MAX) < Tm) {
//         s->y_bin ^= (1 << (rand() % BITS)); // Inverte um bit aleatório de y
//     }
// }
//
// // Função para mutação randômica não uniforme
// void mutate_non_uniform(Solution* s, double Tm, int generation, int max_generations) {
//     if (((double)rand() / RAND_MAX) < Tm * (1 - (double)generation / max_generations)) {
//         s->x_bin ^= (1 << (rand() % BITS)); // Inverte um bit aleatório de x
//     }
//     if (((double)rand() / RAND_MAX) < Tm * (1 - (double)generation / max_generations)) {
//         s->y_bin ^= (1 << (rand() % BITS)); // Inverte um bit aleatório de y
//     }
// }

// Função que executa o algoritmo genético
void run_genetic_algorithm(int population_size, int generations, const char* filename, double Tc, double Tm,
                           void (*crossover_func)(Solution, Solution, Solution*, Solution*, double),
                           void (*mutate_func)(Solution*, double)) {
    Solution* population = (Solution*)malloc(population_size * sizeof(Solution));
    if (population == NULL) {
        fprintf(stderr, "Erro ao alocar memória.\n");
        return;
    }

    srand((unsigned int)time(NULL));

    // Inicializa a população com valores binários aleatórios
    for (int i = 0; i < population_size; i++) {
        population[i].x_bin = rand() % (1 << BITS); // Gera um valor binário aleatório para x
        population[i].y_bin = rand() % (1 << BITS); // Gera um valor binário aleatório para y
        calculate_fitness(&population[i]); // Calcula a aptidão da solução
    }

    FILE* output_file = fopen(filename, "w");
    if (!output_file) {
        fprintf(stderr, "Erro ao abrir o arquivo para escrita.\n");
        free(population);
        return;
    }


    // // Variáveis para rastrear o melhor fitness
    double best_fitness = -INFINITY; // Inicializa com o menor valor possível
    double best_x = 0.0;
    double best_y = 0.0;
    int best_generation = 0;

    // Salva os valores de todos os indivíduos a cada 10 gerações
    FILE* f6_file = fopen("f6_geracoes_bin.txt", "a");
    if (!f6_file) {
        fprintf(stderr, "Erro ao abrir o arquivo f6_geracoes.txt para escrita.\n");
        fclose(output_file);
        free(population);
        return;
    }

    for (int generation = 0; generation < generations; ++generation) {
        Solution* new_population = (Solution*)malloc(population_size * sizeof(Solution));
        if (new_population == NULL) {
            fprintf(stderr, "Erro ao alocar memória para nova população.\n");
            fclose(output_file);
            fclose(f6_file);
            free(population);
            return;
        }

        for (int i = 0; i < population_size; i += 2) {
            Solution parent1 = roulette_selection(population, population_size);
            Solution parent2 = roulette_selection(population, population_size);
            Solution child1, child2;

            crossover_one_point(parent1, parent2, &child1, &child2, Tc);
            mutate(&child1, Tm);
            mutate(&child2, Tm);

            calculate_fitness(&child1);
            calculate_fitness(&child2);

            new_population[i] = child1;
            new_population[i + 1] = child2;
        }

        for (int i = 0; i < population_size; i++) {
            population[i] = new_population[i];
        }

        // Salva os valores de todos os indivíduos a cada 10 gerações
        if (generation % 10 == 0) {
            for (int j = 0; j < population_size; j++) {
                fprintf(f6_file, "%d %f %f %f\n", generation, decode(population[j].x_bin), decode(population[j].y_bin), population[j].fitness_value);
            }
        }

        Solution best_solution = population[0];
        for (int i = 1; i < population_size; i++) {
            if (population[i].fitness_value > best_solution.fitness_value) {
                best_solution = population[i];
            }
        }

        // Atualiza o melhor fitness encontrado
        if (best_solution.fitness_value > best_fitness) {
            best_fitness = best_solution.fitness_value;
            best_x = decode(best_solution.x_bin);
            best_y = decode(best_solution.y_bin);
            best_generation = generation; // Armazena a geração do melhor fitness
        }


        // Saída para o arquivo
        fprintf(output_file, "%d %f %f %f\n",
                generation, best_solution.fitness_value, decode(best_solution.x_bin), decode(best_solution.y_bin));

        // Impressão no console
        printf("Geração %d: Melhor Fitness = %f | x = %f | y = %f\n",
            generation, best_solution.fitness_value, decode(best_solution.x_bin), decode(best_solution.y_bin));


        // Salva a matriz binária a cada 10 gerações
        if (generation % 10 == 0) {
            char x_bin_pop[BITS + 1]; // +1 para o caractere nulo
            char y_bin_pop[BITS + 1]; // +1 para o caractere nulo
            FILE* matrizs = fopen("matrizes.txt", "a"); // Abre o arquivo em modo de anexação
            if (!matrizs) {
                fprintf(stderr, "Erro ao abrir o arquivo para escrita.\n");
                free(new_population);
                continue; // Continua para a próxima geração
            }

            for (int j = 0; j < population_size; j++) {
                int_to_binary_string(population[j].x_bin, x_bin_pop, sizeof(x_bin_pop));
                int_to_binary_string(population[j].y_bin, y_bin_pop, sizeof(y_bin_pop));
                fprintf(matrizs, "%d %s %s\n", generation, x_bin_pop, y_bin_pop); // Salva as representações binárias
            }
            fclose(matrizs); // Fecha o arquivo após a escrita
        }

        // Converte e imprime a representação binária
        char x_bin_str[BITS + 1]; // +1 para o caractere nulo
        char y_bin_str[BITS + 1]; // +1 para o caractere nulo
        int_to_binary_string(best_solution.x_bin, x_bin_str, sizeof(x_bin_str));
        int_to_binary_string(best_solution.y_bin, y_bin_str, sizeof(y_bin_str));
        printf("Binário = %s , %s\n", x_bin_str, y_bin_str);

        free(new_population);

    }

    // // Salva o melhor fitness encontrado e suas coordenadas
    // fprintf(output_file, "%d %f %f %f\n",
    //         best_generation, best_fitness, best_x, best_y);
    //
    printf("Melhor Fitness encontrado: %f em Geração %d | x = %f | y = %f\n",
             best_fitness, best_generation, best_x, best_y);

    fclose(f6_file); // Fecha o arquivo após a escrita
    fclose(output_file);
    free(population);
}

void teste(double x, double y) {
    double z = 0.5 - (pow(sin(sqrt(x * x + y * y)), 2) - 0.5) / pow(1 + 0.001 * (x * x + y * y), 2);
    printf("Resultado para x = %f e y = %f: %f", x, y, z);
}

int main() {
    run_genetic_algorithm(200, 500, "results_bin.txt", 0.85, 0.01, crossover_one_point, mutate);
    return 0; // Retorna 0 para indicar que o programa terminou com sucesso
}
