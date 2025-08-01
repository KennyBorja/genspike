# Instrucciones de Ejecución - GenSpike

## Resumen del Proyecto

GenSpike es una aplicación web Django que implementa un algoritmo genético para simular la evolución del gen Spike del SARS-CoV-2. El proyecto integra datos de múltiples fuentes biológicas (NCBI, UniProt, EBI) para proporcionar análisis evolutivo basado en evidencia experimental.

## Requisitos Previos

### Software Necesario
- Python 3.11 o superior
- pip (gestor de paquetes de Python)
- Navegador web moderno
- Conexión a internet estable

### Recursos del Sistema
- RAM: 4GB mínimo, 8GB recomendado
- Espacio en disco: 1GB libre
- CPU: Procesador moderno (recomendado multi-core)

## Instalación Paso a Paso

### 1. Preparación del Entorno

```bash
# Navegar al directorio del proyecto
cd genspike_project

# (Opcional) Crear entorno virtual
python -m venv genspike_env
source genspike_env/bin/activate  # En Linux/Mac
# genspike_env\Scripts\activate   # En Windows
```

### 2. Instalación de Dependencias

```bash
# Instalar todas las dependencias
pip install -r requirements.txt

# Verificar instalación
python -c "import django; print(django.get_version())"
```

### 3. Configuración de Django

```bash
# Ejecutar migraciones de base de datos
python manage.py migrate

# (Opcional) Crear superusuario para admin
python manage.py createsuperuser
```

## Ejecución de la Aplicación

### Iniciar el Servidor

```bash
# Iniciar servidor de desarrollo
python manage.py runserver 0.0.0.0:8000

# El servidor estará disponible en:
# http://localhost:8000
# http://127.0.0.1:8000
```

### Acceso a la Aplicación

1. Abrir navegador web
2. Navegar a `http://localhost:8000`
3. La interfaz principal se cargará automáticamente

## Uso de la Aplicación

### Flujo de Trabajo Recomendado

1. **Cargar Datos Iniciales**
   - Hacer clic en "Cargar Gen Spike"
   - Hacer clic en "Cargar Proteína Original"
   - Esperar confirmación de carga exitosa

2. **Configurar Parámetros del Algoritmo**
   - Tamaño de Población: 50 (recomendado)
   - Generaciones: 20 (recomendado)
   - Tasa de Mutación: 0.01 (recomendado)

3. **Ejecutar Simulación**
   - Hacer clic en "Inicializar Población"
   - Hacer clic en "Ejecutar Evolución"
   - Monitorear progreso en tiempo real

4. **Analizar Resultados**
   - Revisar gráficos de evolución de fitness
   - Examinar análisis del mejor individuo
   - Explorar mutaciones identificadas

### Parámetros Avanzados

- **Tamaño de Población**: Controla diversidad vs. velocidad
  - Pequeño (20-30): Rápido pero menos diverso
  - Medio (50-80): Balance recomendado
  - Grande (100+): Más diverso pero más lento

- **Número de Generaciones**: Controla convergencia
  - Pocas (5-10): Resultados preliminares
  - Medio (15-25): Análisis estándar
  - Muchas (30+): Análisis exhaustivo

- **Tasa de Mutación**: Controla exploración
  - Baja (0.001-0.005): Conservativa
  - Media (0.01): Recomendada
  - Alta (0.05+): Exploratoria

## APIs Disponibles

### Endpoints Principales

- `GET /api/extract-spike-gene/` - Extraer gen Spike
- `GET /api/get-original-protein/` - Obtener proteína original
- `GET /api/get-functional-variants/` - Obtener variantes funcionales
- `POST /api/initialize-population/` - Inicializar población
- `POST /api/run-evolution/` - Ejecutar algoritmo genético
- `GET /api/analyze-best-individual/` - Analizar mejor individuo
- `GET /api/get-population-status/` - Estado de población
- `POST /api/reset-population/` - Resetear población

### Ejemplo de Uso de API

```python
import requests

# Inicializar población
response = requests.post('http://localhost:8000/api/initialize-population/', 
                        json={'population_size': 50})

# Ejecutar evolución
response = requests.post('http://localhost:8000/api/run-evolution/', 
                        json={
                            'population_size': 50,
                            'generations': 20,
                            'mutation_rate': 0.01
                        })

# Obtener resultados
results = response.json()
```

## Solución de Problemas

### Errores Comunes

1. **Error de Conexión a APIs Externas**
   - Verificar conexión a internet
   - Reintentar después de unos minutos
   - Algunas APIs pueden tener límites de tasa

2. **Rendimiento Lento**
   - Reducir tamaño de población
   - Reducir número de generaciones
   - Cerrar otras aplicaciones

3. **Errores de Memoria**
   - Reducir parámetros del algoritmo
   - Reiniciar el servidor
   - Verificar recursos del sistema

4. **Errores de Instalación**
   - Verificar versión de Python
   - Actualizar pip: `pip install --upgrade pip`
   - Instalar dependencias individualmente

### Logs y Debugging

```bash
# Ejecutar con logs detallados
python manage.py runserver --verbosity=2

# Ver logs en tiempo real
tail -f /var/log/django.log  # Si está configurado
```

## Estructura del Proyecto

```
genspike_project/
├── genetic_algorithm/              # Aplicación principal
│   ├── services/                   # Lógica de negocio
│   │   ├── spike_extractor.py     # Extracción gen Spike
│   │   ├── protein_service.py     # Manejo de proteínas
│   │   ├── variant_service.py     # Variantes funcionales
│   │   └── genetic_algorithm_service.py  # Algoritmo principal
│   ├── templates/                  # Templates HTML
│   ├── static/                     # Archivos estáticos
│   ├── data/                       # Datos de ejemplo
│   └── utils/                      # Utilidades
├── genspike_web/                   # Configuración Django
├── requirements.txt                # Dependencias
├── README.md                       # Documentación completa
├── pipeline_diagram.png            # Diagrama del pipeline
└── genetic_algorithm_pipeline.png  # Diagrama técnico
```

## Personalización y Extensión

### Modificar Parámetros por Defecto

Editar `genetic_algorithm/services/genetic_algorithm_service.py`:

```python
# Cambiar parámetros por defecto
self.population_size = 100  # Cambiar de 50 a 100
self.mutation_rate = 0.005  # Cambiar de 0.01 a 0.005
self.max_generations = 30   # Cambiar de 20 a 30
```

### Agregar Nuevas Métricas de Fitness

1. Modificar `calculate_fitness()` en `GeneticAlgorithmService`
2. Agregar nuevos servicios en `genetic_algorithm/services/`
3. Actualizar la interfaz web según sea necesario

### Integrar Nuevas Fuentes de Datos

1. Crear nuevo servicio en `genetic_algorithm/services/`
2. Implementar métodos de extracción y caché
3. Integrar en el cálculo de fitness

## Contacto y Soporte

Para preguntas técnicas o reportar problemas:
- Revisar la documentación completa en `README.md`
- Verificar la estructura del código en los servicios
- Consultar los diagramas de pipeline incluidos

## Licencia y Uso

Este proyecto está desarrollado con fines educativos y de investigación. Se recomienda citar apropiadamente si se utiliza en publicaciones académicas.

