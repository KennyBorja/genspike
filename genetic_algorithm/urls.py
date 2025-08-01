"""
URLs para la aplicación genetic_algorithm
"""
from django.urls import path
from . import views

app_name = 'genetic_algorithm'

urlpatterns = [
    # Vista principal
    path('', views.index, name='index'),
    
    # APIs para datos básicos
    path('api/extract-spike-gene/', views.extract_spike_gene, name='extract_spike_gene'),
    path('api/get-original-protein/', views.get_original_protein, name='get_original_protein'),
    path('api/get-functional-variants/', views.get_functional_variants, name='get_functional_variants'),
    
    # APIs del algoritmo genético
    path('api/initialize-population/', views.initialize_population, name='initialize_population'),
    path('api/run-evolution/', views.run_evolution, name='run_evolution'),
    path('api/analyze-best-individual/', views.analyze_best_individual, name='analyze_best_individual'),
    
    # APIs de estado
    path('api/get-population-status/', views.get_population_status, name='get_population_status'),
    path('api/reset-population/', views.reset_population, name='reset_population'),
]

