from __future__ import absolute_import

from . import ZmqSubscriber

subscriber = ZmqSubscriber()
subscriber.listen_and_save()