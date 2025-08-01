"""
Servicio para extraer la secuencia del gen Spike del SARS-CoV-2
"""
from Bio import Entrez, SeqIO
import logging

logger = logging.getLogger(__name__)

class SpikeExtractor:
    def __init__(self):
        Entrez.email = "genspike@example.com"
        self.spike_sequence = None
        self.original_record = None
    
    def extract_spike_gene(self):
        """
        Extrae la secuencia del gen Spike del genoma de SARS-CoV-2
        """
        try:
            # Obtener el genoma completo de SARS-CoV-2
            handle = Entrez.efetch(db="nucleotide", id="NC_045512.2", rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            handle.close()
            
            self.original_record = record
            
            # Extraer gen S (Spike) - coordenadas del gen Spike
            start = 21562  # Python usa índice 0, así que restamos 1
            end = 25384
            spike_gene = record.seq[start:end]
            
            self.spike_sequence = spike_gene
            
            logger.info(f"Gen Spike extraído exitosamente. Longitud: {len(spike_gene)} nucleótidos")
            
            return {
                'sequence': str(spike_gene),
                'length': len(spike_gene),
                'start_position': start,
                'end_position': end,
                'first_90': str(spike_gene[:90]),
                'last_90': str(spike_gene[-90:])
            }
            
        except Exception as e:
            logger.error(f"Error al extraer el gen Spike: {str(e)}")
            raise Exception(f"Error al extraer el gen Spike: {str(e)}")
    
    def get_cached_spike_sequence(self):
        """
        Retorna la secuencia del gen Spike si ya fue extraída
        """
        if self.spike_sequence is None:
            return self.extract_spike_gene()
        
        return {
            'sequence': str(self.spike_sequence),
            'length': len(self.spike_sequence),
            'start_position': 21562,
            'end_position': 25384,
            'first_90': str(self.spike_sequence[:90]),
            'last_90': str(self.spike_sequence[-90:])
        }

