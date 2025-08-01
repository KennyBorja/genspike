"""
Vistas para la API del algoritmo genético
"""
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
from django.views.decorators.http import require_http_methods
from django.shortcuts import render
import json
import logging

from .services.genetic_algorithm_service import GeneticAlgorithmService
from .services.spike_extractor import SpikeExtractor
from .services.protein_service import ProteinService
from .services.variant_service import VariantService

logger = logging.getLogger(__name__)

# Instancia global del servicio (simulando estado en memoria)
genetic_service = GeneticAlgorithmService()

def index(request):
    """
    Vista principal con interfaz web
    """
    return render(request, 'genetic_algorithm/index.html')

@csrf_exempt
@require_http_methods(["GET"])
def extract_spike_gene(request):
    """
    API para extraer el gen Spike
    """
    try:
        spike_extractor = SpikeExtractor()
        result = spike_extractor.extract_spike_gene()
        
        return JsonResponse({
            'success': True,
            'data': result
        })
    except Exception as e:
        logger.error(f"Error en extract_spike_gene: {str(e)}")
        return JsonResponse({
            'success': False,
            'error': str(e)
        }, status=500)

@csrf_exempt
@require_http_methods(["GET"])
def get_original_protein(request):
    """
    API para obtener la proteína original
    """
    try:
        protein_service = ProteinService()
        result = protein_service.get_original_protein()
        
        return JsonResponse({
            'success': True,
            'data': result
        })
    except Exception as e:
        logger.error(f"Error en get_original_protein: {str(e)}")
        return JsonResponse({
            'success': False,
            'error': str(e)
        }, status=500)

@csrf_exempt
@require_http_methods(["GET"])
def get_functional_variants(request):
    """
    API para obtener variantes funcionales
    """
    try:
        variant_service = VariantService()
        result = variant_service.get_functional_variants()
        
        # Convertir DataFrame a dict para JSON
        result['dataframe'] = result['dataframe'].to_dict('records')
        
        return JsonResponse({
            'success': True,
            'data': result
        })
    except Exception as e:
        logger.error(f"Error en get_functional_variants: {str(e)}")
        return JsonResponse({
            'success': False,
            'error': str(e)
        }, status=500)

@csrf_exempt
@require_http_methods(["POST"])
def initialize_population(request):
    """
    API para inicializar la población
    """
    try:
        data = json.loads(request.body) if request.body else {}
        population_size = data.get('population_size', 50)
        
        result = genetic_service.initialize_population(population_size)
        
        return JsonResponse({
            'success': True,
            'data': result
        })
    except Exception as e:
        logger.error(f"Error en initialize_population: {str(e)}")
        return JsonResponse({
            'success': False,
            'error': str(e)
        }, status=500)

@csrf_exempt
@require_http_methods(["POST"])
def run_evolution(request):
    """
    API para ejecutar el algoritmo genético
    """
    try:
        data = json.loads(request.body) if request.body else {}
        generations = data.get('generations', 20)
        population_size = data.get('population_size', 50)
        mutation_rate = data.get('mutation_rate', 0.01)
        
        # Configurar parámetros
        genetic_service.max_generations = generations
        genetic_service.population_size = population_size
        genetic_service.mutation_rate = mutation_rate
        
        result = genetic_service.run_evolution(generations)
        
        # Convertir objetos complejos para JSON
        response_data = {
            'best_individual': result['best_individual'],
            'generation_history': result['generation_history'],
            'population_size': len(result['final_population'])
        }
        
        return JsonResponse({
            'success': True,
            'data': response_data
        })
    except Exception as e:
        logger.error(f"Error en run_evolution: {str(e)}")
        return JsonResponse({
            'success': False,
            'error': str(e)
        }, status=500)

@csrf_exempt
@require_http_methods(["GET"])
def analyze_best_individual(request):
    """
    API para analizar el mejor individuo
    """
    try:
        if not genetic_service.population:
            return JsonResponse({
                'success': False,
                'error': 'No hay población disponible. Ejecute primero el algoritmo genético.'
            }, status=400)
        
        # Obtener el mejor individuo
        best_individual = max(genetic_service.population, key=lambda x: x['fitness'])
        
        # Analizar el mejor individuo
        analysis = genetic_service.analyze_best_individual(best_individual)
        
        return JsonResponse({
            'success': True,
            'data': analysis
        })
    except Exception as e:
        logger.error(f"Error en analyze_best_individual: {str(e)}")
        return JsonResponse({
            'success': False,
            'error': str(e)
        }, status=500)

@csrf_exempt
@require_http_methods(["GET"])
def get_population_status(request):
    """
    API para obtener el estado actual de la población
    """
    try:
        if not genetic_service.population:
            return JsonResponse({
                'success': True,
                'data': {
                    'population_exists': False,
                    'population_size': 0,
                    'generations_run': 0
                }
            })
        
        # Calcular estadísticas actuales
        fitness_values = [ind['fitness'] for ind in genetic_service.population]
        
        return JsonResponse({
            'success': True,
            'data': {
                'population_exists': True,
                'population_size': len(genetic_service.population),
                'generations_run': len(genetic_service.generation_history),
                'best_fitness': max(fitness_values) if fitness_values else 0,
                'average_fitness': sum(fitness_values) / len(fitness_values) if fitness_values else 0,
                'worst_fitness': min(fitness_values) if fitness_values else 0
            }
        })
    except Exception as e:
        logger.error(f"Error en get_population_status: {str(e)}")
        return JsonResponse({
            'success': False,
            'error': str(e)
        }, status=500)

@csrf_exempt
@require_http_methods(["POST"])
def reset_population(request):
    """
    API para resetear la población
    """
    try:
        genetic_service.population = []
        genetic_service.generation_history = []
        
        return JsonResponse({
            'success': True,
            'data': {
                'message': 'Población reseteada exitosamente'
            }
        })
    except Exception as e:
        logger.error(f"Error en reset_population: {str(e)}")
        return JsonResponse({
            'success': False,
            'error': str(e)
        }, status=500)

