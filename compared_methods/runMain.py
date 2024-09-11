import sys
import os
import time
import psutil

def measure_performance(func):
    def wrapper(*args, **kwargs):
        process = psutil.Process(os.getpid())
        start_time = time.time()
        start_memory = process.memory_info().rss / 1024 / 1024  # 转换为MB

        result = func(*args, **kwargs)

        end_time = time.time()
        end_memory = process.memory_info().rss / 1024 / 1024  # 转换为MB

        print(f"{func.__name__}:")
        print(f"Runing time: {end_time - start_time:.2f} seconds")
        print(f"Memory usgaes: {end_memory - start_memory:.2f} MB")

        return result
    return wrapper

class Deconvolutions:
    def __init__(self, RNA_file = None, RNA_h5ad = None, RNA_h5Seurat = None, Spatial_file = None, Spatial_h5ad = None, Spatial_h5Seurat = None, celltype_key = None, celltype_file = None, python_path = None, output_path = None):
        """
            @author: wen zhang
            This function integrates spatial and scRNA-seq data to predictes the celltype deconvolution of the spots.
            
            A minimal example usage:
            Assume we have (1) scRNA-seq data file named RNA_h5ad or RNA_h5Seurat
            (2) spatial transcriptomics data file named Spatial_h5ad or Spatial_h5Seurat
            (3) celltype annotataion title in scRNA-seq data file
            
            >>> import Benchmarking.DeconvolutionSpot as DeconvolutionSpot
            >>> test = DeconvolutionSpot.Deconvolutions(RNA_file, RNA_h5ad, RNA_h5Seurat, Spatial_file, Spatial_h5ad, Spatial_h5Seurat, celltype_key, celltype_file, python_path, output_path)
            >>> Methods = ['Cell2location','SpatialDWLS','RCTD','STRIDE','Stereoscope','Tangram','DestVI', 'Seurat', 'SPOTlight', 'DSTG']
            >>> Result = test.Dencon(Methods)
            
            Parameters
            -------
            
            RNA_file : str
            scRNA-seq data count file.
            
            RNA_h5ad : str
            scRNA-seq data file with h5ad format.
            
            RNA_h5Seurat : str
            scRNA-seq data file with h5Seurat format.
            
            Spatial_file : str
            Spatial data count file.
            
            Spatial_h5ad : str
            Spatial data file with h5ad format.
            
            Spatial_h5Seurat : str
            Spatial data file with h5Seurat format.
            
            celltype_key : str
            celltype annotataion title in scRNA-seq data h5ad file or h5Seurat file
            
            celltype_file : str
            celltype annotataion file
            
            python_path : str
            which python path used for Cell2location
            
            output_path : str
            Outfile path
            
            """
        
        self.RNA_file = RNA_file
        self.RNA_h5ad = RNA_h5ad
        self.RNA_h5Seurat = RNA_h5Seurat
        self.Spatial_file = Spatial_file
        self.Spatial_h5ad = Spatial_h5ad
        self.Spatial_h5Seurat = Spatial_h5Seurat
        self.celltype_key = celltype_key
        self.celltype_file = celltype_file
        self.python_path = python_path
        self.output_path = output_path
    
    @measure_performance
    def Dencon(self, need_tools):
        if "Cell2location" in need_tools: #
            RNA_h5ad = self.RNA_h5ad
            Spatial_h5ad = self.Spatial_h5ad
            celltype_key = self.celltype_key
            output_path = self.output_path
            import subprocess
            commands = f'''
                eval "$(conda shell.bash hook)"
                conda activate cell2loc_env
                python tools/Cell2location_pipeline.py {RNA_h5ad} {Spatial_h5ad} {celltype_key} {output_path}
            '''
            subprocess.run(commands, shell=True, executable='/bin/bash')
            #os.system('python tools/Cell2location_pipeline.py ' + RNA_h5ad + ' ' + Spatial_h5ad + ' ' + celltype_key + ' ' + output_path)

        if "SpatialDWLS" in need_tools: #
            RNA_h5Seurat = self.RNA_h5Seurat
            Spatial_h5Seurat = self.Spatial_h5Seurat
            celltype_key = self.celltype_key
            output_path = self.output_path
            python_path = '/root/miniconda3/envs/Benchmarking/bin/python'
            os.system('Rscript tools/SpatialDWLS_pipeline.r ' + RNA_h5Seurat + ' ' + Spatial_h5Seurat + ' ' + celltype_key + ' ' +  python_path + ' ' + output_path)

        if "RCTD" in need_tools: #
            RNA_h5Seurat = self.RNA_h5Seurat
            Spatial_h5Seurat = self.Spatial_h5Seurat
            celltype_key = self.celltype_key
            output_path = self.output_path
            os.system('Rscript tools/RCTD_pipeline.r ' + RNA_h5Seurat + ' ' + Spatial_h5Seurat + ' ' + celltype_key + ' ' + output_path)

        if "STRIDE" in need_tools:
            RNA_file = self.RNA_file
            Spatial_file = self.Spatial_file
            celltype_file = self.celltype_file
            output_path = self.output_path
            os.system('sh tools/STRIDE_pipeline.sh ' + RNA_file + ' ' + Spatial_file + ' ' + celltype_file + ' ' + output_path)

        if "Stereoscope" in need_tools: #
            RNA_h5ad = self.RNA_h5ad
            Spatial_h5ad = self.Spatial_h5ad
            celltype_key = self.celltype_key
            output_path = self.output_path
            os.system('python tools/Stereoscope_pipeline.py ' + RNA_h5ad + ' ' + Spatial_h5ad + ' ' + celltype_key + ' ' + output_path)

        if "Tangram" in need_tools: #
            RNA_h5ad = self.RNA_h5ad
            Spatial_h5ad = self.Spatial_h5ad
            celltype_key = self.celltype_key
            output_path = self.output_path
            os.system('python tools/Tangram_pipeline.py ' + RNA_h5ad + ' ' + Spatial_h5ad + ' ' + celltype_key + ' ' + output_path)

        if "DestVI" in need_tools: #
            RNA_h5ad = self.RNA_h5ad
            Spatial_h5ad = self.Spatial_h5ad
            celltype_key = self.celltype_key
            output_path = self.output_path
            os.system('python tools/DestVI_pipeline.py ' + RNA_h5ad + ' ' + Spatial_h5ad + ' ' + celltype_key + ' ' + output_path)

        if "SpaOTsc" in need_tools: #
            RNA_h5ad = self.RNA_h5ad
            Spatial_h5ad = self.Spatial_h5ad
            celltype_key = self.celltype_key
            output_path = self.output_path
            os.system('python tools/SpaOTsc_pipeline.py ' + RNA_h5ad + ' ' + Spatial_h5ad + ' ' + celltype_key + ' ' + output_path)

        if "novoSpaRc" in need_tools: #
            RNA_h5ad = self.RNA_h5ad
            Spatial_h5ad = self.Spatial_h5ad
            celltype_key = self.celltype_key
            output_path = self.output_path
            os.system('python tools/novoSpaRc_pipeline.py ' + RNA_h5ad + ' ' + Spatial_h5ad + ' ' + celltype_key + ' ' + output_path)

        if "Seurat" in need_tools: #
            RNA_h5Seurat = self.RNA_h5Seurat
            Spatial_h5Seurat = self.Spatial_h5Seurat
            celltype_key = self.celltype_key
            output_path = self.output_path
            os.system('Rscript tools/seurat_pipeline.R ' + RNA_h5Seurat + ' ' + Spatial_h5Seurat + ' ' + celltype_key + ' ' + output_path)

        if "SPOTlight" in need_tools: #
            RNA_h5Seurat = self.RNA_h5Seurat
            Spatial_h5Seurat = self.Spatial_h5Seurat
            celltype_key = self.celltype_key
            output_path = self.output_path
            os.system('Rscript tools/SPOTlight_pipeline.r ' + RNA_h5Seurat + ' ' + Spatial_h5Seurat + ' ' + celltype_key + ' ' + output_path)

        if "CytoSPACE" in need_tools: #
            RNA_h5Seurat = self.RNA_h5Seurat
            Spatial_h5Seurat = self.Spatial_h5Seurat
            celltype_key = self.celltype_key
            output_path = self.output_path
            import subprocess
            commands = f'''
                eval "$(conda shell.bash hook)"
                conda activate cytospace
                Rscript tools/CytoSpace_pipeline.R {RNA_h5Seurat} {Spatial_h5Seurat} {celltype_key} {output_path}
            '''
            subprocess.run(commands, shell=True, executable='/bin/bash')

        if "CARD" in need_tools: #
            RNA_h5Seurat = self.RNA_h5Seurat
            Spatial_h5Seurat = self.Spatial_h5Seurat
            celltype_key = self.celltype_key
            output_path = self.output_path
            os.system('Rscript tools/CARD_pipeline.R ' + RNA_h5Seurat + ' ' + Spatial_h5Seurat + ' ' + celltype_key + ' ' + output_path)

        if "CellTrek" in need_tools: #
            RNA_h5Seurat = self.RNA_h5Seurat
            Spatial_h5Seurat = self.Spatial_h5Seurat
            celltype_key = self.celltype_key
            output_path = self.output_path
            import subprocess
            commands = f'''
                eval "$(conda shell.bash hook)"
                conda activate base
                Rscript tools/CellTrek_pipeline.R  {RNA_h5Seurat} {Spatial_h5Seurat} {celltype_key} {output_path}
            '''
            subprocess.run(commands, shell=True, executable='/bin/bash')
            #os.system('Rscript tools/CellTrek_pipeline.R ' + RNA_h5Seurat + ' ' + Spatial_h5Seurat + ' ' + celltype_key + ' ' + output_path)
        
        if "DSTG" in need_tools: #
            RNA_h5Seurat = self.RNA_h5Seurat
            Spatial_h5Seurat = self.Spatial_h5Seurat
            celltype_key = self.celltype_key
            output_path = self.output_path
            prefix = self.celltype_file
            import subprocess
            commands = f'''
                eval "$(conda shell.bash hook)"
                conda activate DSTG
                sh tools/DSTG_pipeline.sh {RNA_h5Seurat} {Spatial_h5Seurat} {celltype_key} {output_path} {prefix}
            '''
            subprocess.run(commands, shell=True, executable='/bin/bash')
            #os.system('sh tools/DSTG_pipeline.sh ' + RNA_h5Seurat + ' ' + Spatial_h5Seurat + ' ' + celltype_key + ' ' + output_path)

        if "Cell2Spatial" in need_tools: #
            RNA_h5Seurat = self.RNA_h5Seurat
            Spatial_h5Seurat = self.Spatial_h5Seurat
            celltype_key = self.celltype_key
            output_path = self.output_path
            import subprocess
            commands = f'''
                eval "$(conda shell.bash hook)"
                conda activate base
                Rscript tools/Cell2Spatial_pipeline.R  {RNA_h5Seurat} {Spatial_h5Seurat} {celltype_key} {output_path}
            '''
            subprocess.run(commands, shell=True, executable='/bin/bash')

if __name__ == "__main__":
    out_dir = 'results2/Sim.0'
    
    rds = ['data/end.to.end/obj.sc.test.RDS', 'data/end.to.end/stxBrain.sim.0.rds', 'CellType', out_dir]
    h5ad = ['data/allen_cortex.h5ad', 'data/stxBrain.sim.5.h5ad', 'subclass', out_dir]
    h5ad = ['data/end.to.end/obj.sc.test.h5ad', 'data/end.to.end/stxBrain.sim.0.h5ad', 'CellType', out_dir]
    fil = ['allen_cortex.txt', 'stxBrain.sim.0.txt', 'cell_types.txt', out_dir]

    test = Deconvolutions(RNA_h5ad = h5ad[0], Spatial_h5ad = h5ad[1], celltype_key = h5ad[2], output_path = h5ad[3])
    Result = test.Dencon("Cell2location")
    Result = test.Dencon("Stereoscope")
    Result = test.Dencon("Tangram")
    Result = test.Dencon("DestVI")
    Result = test.Dencon("SpaOTsc")
    Result = test.Dencon("novoSpaRc")

    test = Deconvolutions(RNA_h5Seurat = rds[0], Spatial_h5Seurat = rds[1], celltype_key = rds[2], output_path = rds[3], celltype_file = "sim0")
    Result = test.Dencon("SpatialDWLS")
    Result = test.Dencon("RCTD")
    Result = test.Dencon("Seurat")
    Result = test.Dencon("SPOTlight")
    Result = test.Dencon("CytoSPACE")
    Result = test.Dencon("CARD")
    Result = test.Dencon("CellTrek")
    Result = test.Dencon("Cell2Spatial")
    Result = test.Dencon("DSTG")

