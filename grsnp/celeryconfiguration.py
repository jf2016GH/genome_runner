from kombu import Exchange, Queue


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

