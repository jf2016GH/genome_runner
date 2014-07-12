from celery import Celery
from celery import signals
import grsnp.hypergeom4
from celery.bin import Option
from celery.exceptions import Reject, MaxRetriesExceededError
import os


# celery
app = Celery('grsnp')
app.config_from_object('grsnp.celeryconfiguration')
app.user_options['preload'].add(
    Option('-d', '--data_dir', default='',
           help='Set the directory containing the database. Required. Use absolute path. Example: /home/username/db_#.##_#.##.####/.'),
)
app.user_options['preload'].add(
    Option('-r', '--run_files_dir', default='',
           help="Set the directory where the server should save results. Required. Use absolute path. Example: /home/username/run_files/."),
)

sett = {}

@app.task(ignore_result=False)
def run_hypergeom(fois, gfs, bg_path,job_name="",zip_run_files=False,bkg_overlaps_path="",run_annotation=False,run_randomization_test=False,padjust="None",pct_score="",organism="",id=""):
	global sett
	outdir=os.path.join(sett['run_files_dir'],'results',str(id))
	# write out absolute gfs and fois file paths and pass these to the worker	
	for f_path in [fois, gfs]:
		list_f = [x for x in open(f_path).read().split("\n") if x!= ""]
		with open(f_path+'_full','wb') as writer:
			for f in list_f:
				# append the database directory or run files directory based on 
				# the whether the path in the .fois or .gfs file starts with /scd_db
				if f.startswith("/scd_db"):
					writer.write(os.path.join(sett['data_dir'],f.lstrip('/'))+"\n")
				else:
					writer.write(os.path.join(sett['data_dir'],f.lstrip('/'))+"\n")
	bg_path = os.path.join(sett['data_dir'],bg_path)
	grsnp.hypergeom4.run_hypergeom(fois+"_full", gfs+"_full", bg_path,outdir,job_name,zip_run_files,bkg_overlaps_path,sett['data_dir'],run_annotation,run_randomization_test,padjust,pct_score,organism)

# process command line arguments if they exist
@signals.user_preload_options.connect
def cmd_options(options,**kwargs):
	# These settings set the location of the database that the celery worker should use.
	# The data_dir can be in a different location from that of the server's data_dir provided
	# it has the exact same GFs as the server.
	# run_files_dir MUST point to the same database that the server is using.
	global sett
	if options['data_dir'] == '':
		raise Exception('data_dir is a required argument')
	if not os.path.exists(options['data_dir']):
		raise Exception('{} does not exist'.format(options['data_dir']))
	if options['run_files_dir'] == '':
		raise Exception('run_files_dir is a required argument')
	if not os.path.exists(options['run_files_dir']):
		raise Exception('{} does not exist'.format(options['run_files_dir']))
	sett["data_dir"] = options['data_dir'].rstrip("/")
	sett["run_files_dir"] = options["run_files_dir"].rstrip("/")

