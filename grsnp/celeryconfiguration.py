from kombu import Exchange, Queue

db_num = "0"

# redis server config code
# redis_port = 7775
# BROKER_URL = "redis://localhost:{}/".format(redis_port) + db_num
6

# Workers should run as an unprivileged user.
CELERYD_USER="ubuntu"
CELERYD_GROUP="ubuntu"


# RabbitMQ configuration
BROKER_URL = "amqp://guest:guest@localhost:5672//"


CELERY_IMPORTS = ("grsnp")

# info about the two settings below can be read at
# http://docs.celeryproject.org/en/latest/userguide/optimizing.html#prefork-pool-prefetch-settings
CELERYD_PREFETCH_MULTIPLIER = 1 
CELERY_TASK_PUBLISH_RETRY = False


#CELERY_SEND_TASK_SENT_EVENT = True


## Celery queue settings
#CELERY_ROUTES = ({'grsnp.worker_hypergeom4.run_hypergeom': {'queue': 'short_runs'}},
#					{'grsnp.worker_hypergeom4.run_hypergeom': {'queue': 'long_runs'}}
#				)

#default_exchange = Exchange('default', type='direct')

#CELERY_QUEUES = (
#    Queue('default', default_exchange, routing_key='default'),
#    Queue('long_runs', default_exchange, routing_key='grsnp'),
#    Queue('short_runs', default_exchange, routing_key='grsnp')
#)
#CELERY_DEFAULT_QUEUE = 'default'
#CELERY_DEFAULT_EXCHANGE = 'default'
#CELERY_DEFAULT_ROUTING_KEY = 'default'
