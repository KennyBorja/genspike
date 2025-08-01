"""
Servicio principal del algoritmo genético para evolución del gen Spike
"""
import random
import copy
import logging
from typing import List, Dict, Any
from .spike_extractor import SpikeExtractor
from .protein_service import ProteinService
from .variant_service import VariantService
from .mutation_analysis_service import MutationAnalysisService

logger = logging.getLogger(__name__)

class GeneticAlgorithmService:
    def __init__(self):
        self.spike_extractor = SpikeExtractor()
        self.protein_service = ProteinService()
        self.variant_service = VariantService()
        self.mutation_analysis_service = MutationAnalysisService()
        
        self.original_sequence = None
        self.original_protein = None
        self.population = []
        self.generation_history = []
        
        # Parámetros del algoritmo genético
        self.population_size = 50
        self.mutation_rate = 0.01
        self.max_generations = 20
        
    def initialize_population(self, population_size=None):
        """
        Inicializa la población con mutaciones aleatorias
        """
        if population_size:
            self.population_size = population_size
            
        try:
            # Obtener la secuencia original del gen Spike
            spike_data = self.spike_extractor.get_cached_spike_sequence()
            self.original_sequence = spike_data['sequence']
            
            # Obtener la proteína original
            protein_data = self.protein_service.get_cached_protein()
            self.original_protein = protein_data['sequence']
            
            # Generar población inicial
            self.population = []
            
            for i in range(self.population_size):
                # Crear individuo basado en la secuencia original con mutaciones aleatorias
                individual = self.create_individual_with_mutations(self.original_sequence)
                self.population.append(individual)
            
            logger.info(f"Población inicial creada con {len(self.population)} individuos")
            
            return {
                'population_size': len(self.population),
                'original_sequence_length': len(self.original_sequence),
                'original_protein_length': len(self.original_protein)
            }
            
        except Exception as e:
            logger.error(f"Error al inicializar la población: {str(e)}")
            raise Exception(f"Error al inicializar la población: {str(e)}")
    
    def create_individual_with_mutations(self, original_sequence):
        """
        Crea un individuo con mutaciones aleatorias
        """
        sequence = list(original_sequence)
        nucleotides = ['A', 'T', 'G', 'C']
        
        # Aplicar mutaciones aleatorias (pequeñas para mantener diversidad inicial)
        num_mutations = random.randint(1, 5)  # Entre 1 y 5 mutaciones iniciales
        
        for _ in range(num_mutations):
            position = random.randint(0, len(sequence) - 1)
            original_nucleotide = sequence[position]
            new_nucleotide = random.choice([n for n in nucleotides if n != original_nucleotide])
            sequence[position] = new_nucleotide
        
        individual = {
            'sequence': ''.join(sequence),
            'fitness': 0.0,
            'nucleotide_fitness': 0.0,
            'functional_fitness': 0.0,
            'protein_sequence': '',
            'mutations': [],
            'generation': 0
        }
        
        return individual
    
    def calculate_fitness(self, individual):
        """
        Calcula el fitness total de un individuo
        """
        try:
            # 1. Calcular fitness nucleotídico (penalizar cambios en ADN)
            nucleotide_fitness = self.calculate_nucleotide_fitness(individual['sequence'])
            
            # 2. Convertir a proteína
            protein_sequence = self.protein_service.translate_nucleotides_to_protein(individual['sequence'])
            individual['protein_sequence'] = protein_sequence
            
            # 3. Comparar con proteína original
            protein_changes = self.protein_service.compare_proteins(self.original_protein, protein_sequence)
            individual['mutations'] = protein_changes
            
            # 4. Calcular fitness funcional
            functional_fitness = self.variant_service.calculate_functional_fitness(protein_changes)
            
            # 5. Combinar ambos fitness
            total_fitness = functional_fitness * nucleotide_fitness
            
            individual['nucleotide_fitness'] = nucleotide_fitness
            individual['functional_fitness'] = functional_fitness
            individual['fitness'] = total_fitness
            
            return individual
            
        except Exception as e:
            logger.error(f"Error al calcular fitness: {str(e)}")
            individual['fitness'] = 0.01  # Fitness mínimo en caso de error
            return individual
    
    def calculate_nucleotide_fitness(self, sequence):
        """
        Calcula el fitness nucleotídico penalizando cambios en el ADN
        """
        if not self.original_sequence:
            return 1.0
        
        differences = sum(1 for a, b in zip(self.original_sequence, sequence) if a != b)
        similarity = 1.0 - (differences / len(self.original_sequence))
        
        # Aplicar una función exponencial para penalizar más los cambios
        fitness = similarity ** 2
        
        return max(fitness, 0.01)  # Mínimo fitness de 0.01
    
    def select_parents(self):
        """
        Selección de padres proporcional al fitness
        """
        # Calcular fitness de toda la población
        for individual in self.population:
            self.calculate_fitness(individual)
        
        # Obtener pesos (fitness) para la selección
        weights = [individual['fitness'] for individual in self.population]
        
        # Seleccionar padres
        parents = random.choices(self.population, weights=weights, k=self.population_size)
        
        return parents
    
    def mutate_individual(self, individual, generation):
        """
        Aplica mutación a un individuo
        """
        sequence = list(individual['sequence'])
        nucleotides = ['A', 'T', 'G', 'C']
        
        # Aplicar mutaciones basadas en la tasa de mutación
        for i in range(len(sequence)):
            if random.random() < self.mutation_rate:
                original_nucleotide = sequence[i]
                new_nucleotide = random.choice([n for n in nucleotides if n != original_nucleotide])
                sequence[i] = new_nucleotide
        
        # Crear nuevo individuo
        new_individual = {
            'sequence': ''.join(sequence),
            'fitness': 0.0,
            'nucleotide_fitness': 0.0,
            'functional_fitness': 0.0,
            'protein_sequence': '',
            'mutations': [],
            'generation': generation
        }
        
        return new_individual
    
    def run_evolution(self, generations=None):
        """
        Ejecuta el algoritmo genético completo
        """
        if generations:
            self.max_generations = generations
        
        try:
            # Inicializar población si no existe
            if not self.population:
                self.initialize_population()
            
            self.generation_history = []
            
            for generation in range(self.max_generations):
                logger.info(f"Ejecutando generación {generation + 1}/{self.max_generations}")
                
                # Seleccionar padres
                parents = self.select_parents()
                
                # Crear nueva generación
                new_population = []
                for parent in parents:
                    child = self.mutate_individual(parent, generation + 1)
                    new_population.append(child)
                
                self.population = new_population
                
                # Calcular fitness de la nueva población
                for individual in self.population:
                    self.calculate_fitness(individual)
                
                # Guardar estadísticas de la generación
                generation_stats = self.get_generation_statistics(generation + 1)
                self.generation_history.append(generation_stats)
                
                logger.info(f"Generación {generation + 1}: Mejor fitness = {generation_stats['best_fitness']:.4f}")
            
            # Seleccionar el mejor individuo
            best_individual = max(self.population, key=lambda x: x['fitness'])
            
            return {
                'best_individual': best_individual,
                'generation_history': self.generation_history,
                'final_population': self.population
            }
            
        except Exception as e:
            logger.error(f"Error durante la evolución: {str(e)}")
            raise Exception(f"Error durante la evolución: {str(e)}")
    
    def get_generation_statistics(self, generation):
        """
        Obtiene estadísticas de una generación
        """
        fitness_values = [ind['fitness'] for ind in self.population]
        nucleotide_fitness_values = [ind['nucleotide_fitness'] for ind in self.population]
        functional_fitness_values = [ind['functional_fitness'] for ind in self.population]
        
        return {
            'generation': generation,
            'best_fitness': max(fitness_values),
            'average_fitness': sum(fitness_values) / len(fitness_values),
            'worst_fitness': min(fitness_values),
            'best_nucleotide_fitness': max(nucleotide_fitness_values),
            'average_nucleotide_fitness': sum(nucleotide_fitness_values) / len(nucleotide_fitness_values),
            'best_functional_fitness': max(functional_fitness_values),
            'average_functional_fitness': sum(functional_fitness_values) / len(functional_fitness_values),
            'population_size': len(self.population)
        }
    
    def analyze_best_individual(self, best_individual):
        """
        Analiza el mejor individuo al final de la evolución
        """
        try:
            analysis = {
                'sequence_comparison': {
                    'original_length': len(self.original_sequence),
                    'final_length': len(best_individual['sequence']),
                    'nucleotide_changes': sum(1 for a, b in zip(self.original_sequence, best_individual['sequence']) if a != b)
                },
                'protein_comparison': {
                    'original_length': len(self.original_protein),
                    'final_length': len(best_individual['protein_sequence']),
                    'amino_acid_changes': len(best_individual['mutations'])
                },
                'fitness_breakdown': {
                    'total_fitness': best_individual['fitness'],
                    'nucleotide_fitness': best_individual['nucleotide_fitness'],
                    'functional_fitness': best_individual['functional_fitness']
                },
                'mutations': best_individual['mutations'],
                'generation': best_individual['generation']
            }

            # Realizar análisis de mutaciones adicionales
            evolved_sequences = [best_individual['sequence']]
            mutation_analysis_results = self.mutation_analysis_service.analyze_mutations(
                self.original_sequence,
                evolved_sequences,
                self.max_generations # Usar el número de generaciones del algoritmo genético
            )
            if mutation_analysis_results:
                analysis['mutation_analysis'] = mutation_analysis_results

            # Obtener información detallada de las mutaciones
            detailed_mutations = []
            for mutation in best_individual['mutations']:
                if mutation['position'] != 'length':
                    mutation_info = self.variant_service.get_mutation_info(
                        mutation['position'],
                        mutation['original'],
                        mutation['mutated']
                    )
                    detailed_mutations.append(mutation_info)
            analysis['detailed_mutations'] = detailed_mutations
            
            return analysis
            
        except Exception as e:
            logger.error(f"Error al analizar el mejor individuo: {str(e)}")
            raise Exception(f"Error al analizar el mejor individuo: {str(e)}")

