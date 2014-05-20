from celery import Celery
import grsnp.hypergeom4



# celery
app = Celery('grsnp')
app.config_from_object('grsnp.celeryconfiguration')

@app.task(ignore_result=True)
def run_hypergeom(fois, gfs, bg_path,outdir,job_name="",zip_run_files=False,bkg_overlaps_path="",gr_data_dir = "" ,run_annotation=True,run_randomization_test=False,padjust="None",pct_score=""):	
	grsnp.hypergeom4.run_hypergeom(fois, gfs, bg_path,outdir,job_name,zip_run_files,bkg_overlaps_path,gr_data_dir,run_annotation,run_randomization_test,padjust,pct_score)