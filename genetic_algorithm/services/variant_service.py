"""
Servicio para obtener variantes funcionales desde EBI
"""
import requests
import pandas as pd
import logging

logger = logging.getLogger(__name__)

class VariantService:
    def __init__(self):
        self.accession = 'P0DTC2'  # UniProt de la proteína Spike de SARS-CoV-2
        self.variants_data = None
    
    def get_functional_variants(self):
        """
        Obtiene las variantes funcionales desde EBI
        """
        try:
            url = f'https://www.ebi.ac.uk/proteins/api/proteins/{self.accession}'
            resp = requests.get(url, headers={'Accept': 'application/json'})
            
            if resp.status_code == 200:
                seq_data = resp.json()
                full_protein_sequence = seq_data.get('sequence', {}).get('sequence', None)
                variants = seq_data.get('features', [])
                data = []
                for feat in variants:
                    if feat.get('type') == 'VARIANT':
                        begin = feat.get('begin')
                        try:
                            begin_pos = int(begin) if begin is not None else None
                        except (ValueError, TypeError):
                            begin_pos = None
                        

                        original_aaa = full_protein_sequence[begin_pos - 1]

                        
                        data.append({
                            'mutation': f"{original_aaa}{begin}{feat.get('alternativeSequence')}",
                            'position': begin_pos,
                            'original_aa': original_aaa,
                            'alternative_aa': feat.get('alternativeSequence'),
                            'description': feat.get('description', ''),
                            'evidence': feat.get('evidence', 'N/A')
                        })
                
                self.variants_data = pd.DataFrame(data)
                logger.info(f"Variantes funcionales obtenidas: {len(data)} variantes encontradas")
                
                return {
                    'variants': data,
                    'total_variants': len(data),
                    'dataframe': self.variants_data
                }
            else:
                raise Exception(f"Error al obtener variantes: HTTP {resp.status_code}")
                
        except Exception as e:
            logger.error(f"Error al obtener variantes funcionales: {str(e)}")
            raise Exception(f"Error al obtener variantes funcionales: {str(e)}")
    
    def get_cached_variants(self):
        """
        Retorna las variantes si ya fueron obtenidas
        """
        if self.variants_data is None:
            return self.get_functional_variants()
        
        return {
            'variants': self.variants_data.to_dict('records'),
            'total_variants': len(self.variants_data),
            'dataframe': self.variants_data
        }
    
    def calculate_functional_fitness(self, protein_changes):
        """
        Calcula el fitness funcional basado en las variantes conocidas.
        Penaliza de forma lineal en lugar de multiplicativa.
        """
        if self.variants_data is None:
            self.get_functional_variants()

        base_fitness = 1.0
        total_penalty = 0.0

        penalty_known_variant = 0.005     
        penalty_unknown = 0.01  

        for change in protein_changes:
            position = int(change['position'])
            original = change['original'].upper()
            mutated = change['mutated'].upper()

            variant_match = self.variants_data[
                (self.variants_data['position'] == position) &
                (self.variants_data['original_aa'].str.upper() == original) &
                (self.variants_data['alternative_aa'].str.upper() == mutated)
            ]

            if not variant_match.empty:
                total_penalty += penalty_known_variant
            else:
                total_penalty += penalty_unknown

        fitness_score = base_fitness - total_penalty

        return max(fitness_score, 0.01)

    
    def get_mutation_info(self, position, original_aa, mutated_aa):
        """
        Obtiene información sobre una mutación específica
        """
        if self.variants_data is None:
            self.get_functional_variants()


        variant_match = self.variants_data[
            (self.variants_data['position'] == int(position)) &
            (self.variants_data['original_aa'].str.upper() == original_aa.upper()) &
            (self.variants_data['alternative_aa'].str.upper() == mutated_aa.upper())
        ]


        if not variant_match.empty:
            return variant_match.iloc[0].to_dict()
        else:
           
            return {
                'mutation': f"{original_aa}{position}{mutated_aa}",
                'position': position,
                'original_aa': original_aa,
                'alternative_aa': mutated_aa,
                'description': 'no',
                'evidence': 'N/A'
            }

