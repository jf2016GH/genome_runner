BROKER_TRANSPORT = "redis"

redis_port = 7775 # The port number that redis should use
BROKER_HOST = "localhost"  # Maps to redis host.
BROKER_PORT = redis_port         # Maps to redis port.
BROKER_VHOST = "0"         # Maps to database number.

# Workers should run as an unprivileged user.
CELERYD_USER="ubuntu"
CELERYD_GROUP="ubuntu"
CELERY_RESULT_BACKEND = 'redis://localhost:{}/0'.format(redis_port)

CELERY_IMPORTS = ("grsnp")

#CELERY_SEND_TASK_SENT_EVENT = True

#CELERY_ROUTES = ({'grsnp.worker_hypergeom4.run_hypergeom': {'queue': 'short_runs'}},
#					{'grsnp.worker_hypergeom4.run_hypergeom': {'queue': 'long_runs'}}
#				)
