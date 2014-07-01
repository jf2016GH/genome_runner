BROKER_TRANSPORT = "redis"

db_num = "0"         # Maps to database number.
BROKER_URL = "redis://localhost:7775/" + db_num

# Workers should run as an unprivileged user.
CELERYD_USER="ubuntu"
CELERYD_GROUP="ubuntu"
CELERY_RESULT_BACKEND = BROKER_URL

CELERY_IMPORTS = ("grsnp")

# info about the two settings bellow can be read at
# http://docs.celeryproject.org/en/latest/userguide/optimizing.html#prefork-pool-prefetch-settings
CELERYD_PREFETCH_MULTIPLIER = 1

#CELERY_SEND_TASK_SENT_EVENT = True

#CELERY_ROUTES = ({'grsnp.worker_hypergeom4.run_hypergeom': {'queue': 'short_runs'}},
#					{'grsnp.worker_hypergeom4.run_hypergeom': {'queue': 'long_runs'}}
#				)
