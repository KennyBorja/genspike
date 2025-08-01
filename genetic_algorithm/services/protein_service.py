"""
Servicio para obtener la secuencia de proteína Spike desde UniProt
"""
import requests
from Bio import SeqIO
from io import StringIO
import logging

logger = logging.getLogger(__name__)

class ProteinService:
    def __init__(self):
        self.original_protein = None
        self.accession = "P0DTC2"  # UniProt ID para la proteína Spike de SARS-CoV-2
    
    def get_original_protein(self):
        """
        Obtiene la secuencia de proteína original del Spike desde UniProt
        """
        try:
            url = f"https://www.uniprot.org/uniprotkb/{self.accession}.fasta"
            response = requests.get(url)
            
            if response.status_code == 200:
                fasta_str = response.text
                record = SeqIO.read(StringIO(fasta_str), "fasta")
                
                self.original_protein = record.seq
                
                logger.info(f"Proteína Spike obtenida exitosamente. Longitud: {len(record.seq)} aminoácidos")
                
                return {
                    'id': record.id,
                    'sequence': str(record.seq),
                    'length': len(record.seq),
                    'description': record.description
                }
            else:
                raise Exception(f"Error al obtener la proteína: HTTP {response.status_code}")
                
        except Exception as e:
            logger.error(f"Error al obtener la proteína original: {str(e)}")
            raise Exception(f"Error al obtener la proteína original: {str(e)}")
    
    def get_cached_protein(self):
        """
        Retorna la proteína original si ya fue obtenida
        """
        if self.original_protein is None:
            return self.get_original_protein()
        
        return {
            'id': self.accession,
            'sequence': str(self.original_protein),
            'length': len(self.original_protein),
            'description': 'Spike glycoprotein'
        }
    
    def translate_nucleotides_to_protein(self, nucleotide_sequence):
        """
        Traduce una secuencia de nucleótidos a proteína
        """
        try:
            from Bio.Seq import Seq
            
            # Crear objeto Seq y traducir
            seq_obj = Seq(nucleotide_sequence)
            protein_seq = seq_obj.translate()
            
            # Remover el codón de parada si existe
            if protein_seq.endswith('*'):
                protein_seq = protein_seq[:-1]
            
            return str(protein_seq)
            
        except Exception as e:
            logger.error(f"Error al traducir nucleótidos a proteína: {str(e)}")
            raise Exception(f"Error al traducir nucleótidos a proteína: {str(e)}")
    
    def compare_proteins(self, original_protein, mutated_protein):
        """
        Compara dos secuencias de proteínas y retorna las diferencias
        """
        changes = []
        
        min_length = min(len(original_protein), len(mutated_protein))
        
        for i in range(min_length):
            if original_protein[i] != mutated_protein[i]:
                changes.append({
                    'position': i + 1,  # Posición 1-indexada
                    'original': original_protein[i],
                    'mutated': mutated_protein[i]
                })
        
        # Si las secuencias tienen diferentes longitudes
        if len(original_protein) != len(mutated_protein):
            changes.append({
                'position': 'length',
                'original': len(original_protein),
                'mutated': len(mutated_protein)
            })
        
        return changes

