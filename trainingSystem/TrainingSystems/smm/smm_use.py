import subprocess
import tempfile

from smm import SMMTrainInput, TrainingDataType, SequenceDataType, DataPointType, SMMTrainOutput, SMMPredictor, SMMPredictInput, SMMPredictOutput, PredictType, parse

class SMMUse:
    #measurements is a list of the form [sequence, measurement]
    def __init__(self, measurements, alphabet, sequence_length):
        seq_data = SequenceDataType([DataPointType(Sequence = sequence, Measured = measurement) for sequence, measurement in measurements])
        training = TrainingDataType(Alphabet = ''.join(alphabet), SequenceLength = sequence_length, SequenceData = seq_data)
        with tempfile.NamedTemporaryFile() as output:
            smm_train = SMMTrainInput(OutputFile=output.name, TrainingData=training)
            with tempfile.NamedTemporaryFile() as temp:
                with open(temp.name, 'w') as f:
                    smm_train.export(f, 0)
                subprocess.run(['smm', temp.name])
                with open(output.name, 'r') as f:
                    self.predictor = parse(f).get_SMMPredictor()
    """Returns a dictionary mapping a sequence to its predicted measurement"""
    def predict(self, sequences):
        with tempfile.NamedTemporaryFile() as output:
            predict_input = SMMPredictInput(OutputFile = output.name, SMMPredictor = self.predictor,  Predict= [PredictType(sequence) for sequence in sequences])
            with tempfile.NamedTemporaryFile() as prediction_input_file:
                with open(prediction_input_file.name, 'w') as f:
                    predict_input.export(f, 0)
                subprocess.run(['smm', prediction_input_file.name])
                a = output.read()
                with open(output.name, 'r') as f:
                    predict_output = parse(f).get_Predict()
                    return {pOutput.get_Sequence(): float(pOutput.get_Predictions()[0]) for pOutput in predict_output}
"""
measurements = [('DQNP', 2.81), ('KVVK', 0.3), ('YHEL', 1.54), ('HLSV', 3.0), ('ALAV', 2.38), ('ALAV', 1.04), ('AGAV', 2.16), ('YSFL', 1.52), ('AFAF', 0.0), ('ARIA', 0.08), ('IVKK', 1.74), ('KYFY', -0.7), ('RRFF', -1.92), ('LKFI', 0.46), ('LLAI', 1.35), ('VAIE', 2.54), ('RLSI', 3.0), ('GFYV', 1.6), ('AEMY', -0.7), ('QVPY', 0.46), ('ALAV', 1.63), ('KYTL', 0.72), ('FLTL', 2.59), ('FISV', 2.11), ('SRFK', 0.58), ('AVAK', 1.38), ('DLMV', 2.2), ('ALAK', 1.85), ('ADAV', 2.2), ('GLSL', 1.4), ('KKVS', 1.48), ('YQKK', 2.0), ('KRFF', -1.52), ('WLSV', 3.0), ('ALMI', 0.0), ('AGDW', 2.15), ('DVFK', 1.6), ('ARAA', 1.0), ('GQYK', 0.81), ('ALAV', 2.21), ('KIKY', 0.88), ('ARIA', 0.32), ('SRYK', -0.08), ('ALAQ', 2.27), ('LETR', 0.68), ('TKPI', 2.94), ('ALSL', 0.87), ('SAWK', 0.15), ('AAAY', -0.4), ('ANAV', 2.12), ('AYAL', 0.66), ('ALAG', 2.84), ('KPRD', 3.23), ('FADD', 2.26), ('TFGL', 0.7), ('SKSE', 3.07), ('APAM', 1.93), ('ISVF', 1.3), ('TYLL', 0.24), ('AAAY', 0.4), ('ARAA', 1.0), ('ASII', 1.18), ('KYAQ', 2.97), ('SILF', 2.0), ('GRPL', 2.3), ('RYLL', 1.15), ('RRSR', 0.48), ('LLGI', 0.72), ('ATAK', 2.08), ('GLYI', 2.02), ('GEVT', 2.7), ('HEES', 2.41), ('LTRQ', 2.11), ('EIEL', 3.0), ('VMNM', 1.64), ('KYPL', 0.65), ('AIAY', -0.3), ('AAAY', 0.18), ('KWDH', 1.15), ('TRQG', 1.88), ('AIRL', -0.03), ('KIWL', 0.08), ('ALAV', 1.51), ('ATEY', 0.51), ('GEFV', 2.39), ('YLLI', 2.08), ('AAAY', 0.44), ('HEAD', 1.98), ('ALKV', 2.08), ('ALAV', 0.67), ('AVAY', -0.6), ('RRYA', 0.18), ('AKAF', 0.26), ('AAAY', -0.52), ('ILRK', 2.88), ('LLSL', 1.05), ('YFVK', 1.6), ('RAAY', 0.47), ('ARAK', 0.28), ('FLCL', 1.28), ('GLYI', 2.96), ('AYAY', 0.18), ('CLFL', 1.56), ('AAAY', -0.06), ('KYFY', -0.12), ('TITK', 2.67), ('KIKL', 0.93), ('REPM', 2.3), ('QPEY', 2.67), ('YMLR', 2.46), ('KVHY', -0.64), ('ALAT', 2.82), ('IMPI', 3.0), ('AAAY', -0.68), ('ALFL', 0.54), ('RWGQ', 0.7), ('ARAY', -0.38), ('ALFV', 0.43), ('FKDA', 1.81), ('AAEY', 0.43), ('SVVI', 0.54), ('KLHI', 1.78), ('LLDL', 2.31), ('RAYL', -0.2), ('ALAV', 2.21), ('FLWV', 1.46), ('RGQH', 2.18), ('PPEG', 3.25), ('YFDF', 1.38), ('LLFV', 0.9), ('RNVT', 0.82), ('PHIK', 3.07), ('ALAV', 1.41), ('SRVK', 0.7), ('SRWK', -0.47), ('RSQA', 1.93), ('ARLK', 0.65), ('PDVF', 2.9), ('ALAV', 1.52), ('GRQV', 2.7), ('LTLS', 1.7), ('VDRR', 1.89), ('GLFY', -0.05), ('VAFE', 3.12), ('IRKY', -1.0), ('ALAE', 3.27), ('SSKW', 1.18), ('VRRR', 1.32), ('SRYR', -0.62), ('VSSS', 0.48), ('KIKS', 2.38), ('APAF', 1.8), ('RRLL', -0.22), ('AIRL', -0.05), ('GAVF', 0.6), ('IEQK', 1.7), ('RLRH', 2.11), ('EIAL', 2.12), ('AFYA', -0.35), ('ALAV', 1.82), ('RRYL', -0.85), ('AVVT', 1.28), ('SLYL', 0.36), ('CLCV', 2.98), ('ILKV', 2.38), ('DPAD', 2.72), ('AAAY', -0.7), ('DMAR', 1.43), ('ALEV', 2.4), ('QYDF', 1.41), ('PLEL', 2.43), ('QLFI', 3.12), ('WILV', 2.59), ('GRIK', -0.3), ('SRDW', 0.54), ('SLYV', 2.53), ('YDSG', 3.07), ('KAGW', 2.0), ('AQWG', -0.18), ('LLAI', 2.48), ('AADY', -0.1), ('GMRN', 2.73), ('IALR', 0.98), ('YMNV', 2.81), ('AQYL', -0.7), ('QVDL', 0.74), ('HLLL', 2.66), ('ALAL', 0.56), ('GKSL', 1.75), ('EEFA', 1.54), ('YGKV', 2.44), ('RRYL', -0.46), ('KTQF', 0.34), ('LIVA', 2.3), ('PDIS', 3.06), ('YVFD', 3.23), ('EIER', 2.95), ('GRER', 0.94), ('GRAK', 0.08), ('VLGL', 1.86), ('ALMV', 1.4), ('TLWV', 1.4), ('SDII', 1.74), ('ASVP', 1.0), ('YAAY', 1.44), ('KSII', 1.56), ('AAAY', 0.01), ('AYAL', 0.81), ('SWTV', 1.18), ('NQKY', -0.02), ('AFAV', 1.2), ('ALAK', 1.85), ('AYAK', 1.16), ('GIFI', 3.0), ('LSKF', -0.04), ('TNFM', 1.08), ('ALSV', 2.1), ('AARY', -0.21), ('VYSV', 0.81), ('RPQT', 2.12), ('LYNY', 1.38), ('ARAL', 0.28), ('SRFK', 0.34), ('IISV', 1.89), ('AYAK', 1.85), ('ALAA', 2.06), ('GAIK', 1.15), ('AYAF', 0.36), ('FYET', 2.94), ('ALAV', 2.0), ('ERPL', 0.78), ('ALAY', 0.0), ('IVTK', 2.38), ('FVYF', 2.26), ('LYQI', -0.64), ('MVKV', 2.09), ('KLIL', 2.06), ('DSDV', 2.48), ('AADY', 1.21), ('RYRL', 0.74), ('SSKS', 3.11), ('ELDQ', 2.75), ('LVVV', 1.34), ('AAAY', -0.4), ('INLT', 2.97), ('GMAI', 3.0), ('GRIA', 0.97), ('KFKF', 0.4), ('KRYY', -1.59), ('AAAF', 0.0), ('AAAY', 0.02), ('LLQL', 1.39)]
use = SMMUse(measurements, 'ACDEFGHIKLMNPQRSTVWY', 4)
print(use.predict(['AAAA', 'CDEA']))
"""        
