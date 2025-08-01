import numpy as np
from Bio.Seq import Seq
from Bio.Data import CodonTable
from collections import defaultdict
import logging

logger = logging.getLogger(__name__)

class MutationAnalysisService:
    def __init__(self):
        pass

    def calcular_tasa_mutacion(self, ancestral, evolucionadas, generaciones):
        mutaciones_totales = 0
        for sec in evolucionadas:
            mutaciones = sum(1 for a, b in zip(ancestral, sec) if a != b)
            mutaciones_totales += mutaciones
        tasa = mutaciones_totales / (len(ancestral) * len(evolucionadas) / generaciones)
        return tasa

    def espectro_mutaciones(self, ancestral, evolucionadas):
        espectro = defaultdict(int)
        for sec in evolucionadas:
            for i, (a, b) in enumerate(zip(ancestral, sec)):
                if a != b:
                    cambio = f"{a}->{b}"
                    espectro[cambio] += 1
        return dict(espectro)

    def calcular_pN_pS(self, ancestral, evolucionadas):
        standard_table = CodonTable.standard_dna_table
        pN, pS = 0, 0

        for sec in evolucionadas:
            for i in range(0, len(ancestral), 3):
                codon_anc = ancestral[i:i + 3]
                codon_mut = sec[i:i + 3]
                if len(codon_anc) == 3 and len(codon_mut) == 3:
                    try:
                        aa_anc = standard_table.forward_table[str(codon_anc)]
                        aa_mut = standard_table.forward_table[str(codon_mut)]
                        if aa_anc != aa_mut:
                            pN += 1
                        else:
                            pS += 1
                    except KeyError:  # Ignorar codones stop
                        pass
        return pN / pS if pS > 0 else float("inf")

    def mutaciones_recurrentes(self, ancestral, evolucionadas, umbral=2):
        posiciones = defaultdict(int)
        for sec in evolucionadas:
            for i, (a, b) in enumerate(zip(ancestral, sec)):
                if a != b:
                    posiciones[i] += 1
        return {pos: count for pos, count in posiciones.items() if count >= umbral}

    def analyze_mutations(self, ancestral_sequence, evolved_sequences, generations):
        try:
            tasa = self.calcular_tasa_mutacion(ancestral_sequence, evolved_sequences, generations)
            espectro = self.espectro_mutaciones(ancestral_sequence, evolved_sequences)
            pN_pS = self.calcular_pN_pS(ancestral_sequence, evolved_sequences)
            recurrentes = self.mutaciones_recurrentes(ancestral_sequence, evolved_sequences)

            return {
                "tasa_mutacion": tasa,
                "espectro_mutaciones": espectro,
                "pN_pS_ratio": pN_pS,
                "mutaciones_recurrentes": recurrentes
            }
        except Exception as e:
            logger.error(f"Error durante el análisis de mutaciones: {str(e)}")
            return None


