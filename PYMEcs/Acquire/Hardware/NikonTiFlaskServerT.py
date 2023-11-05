from flask import Flask, request, jsonify, make_response
from http import HTTPStatus
import logging

def create_app(mode='simulated'):
    # create and configure the app
    app = Flask(__name__)
    
    if mode == 'simulated':
        from PYMEcs.Acquire.Hardware.LPthreadedSimpleSim import LPThread
    else:
        from PYMEcs.Acquire.Hardware.LPthreadedSimple import LPThread

    # so that we can see which thread processes the various messages, commands
    logging.basicConfig(level=logging.DEBUG,
                        format='(%(threadName)-9s) %(message)s',)

    lpt = LPThread(name='NikonTiThread')
    lpt.start()

    @app.get("/names")
    def get_names():
        status, names = lpt.run_command('GetNames')
        return jsonify(names)

    @app.get("/port")
    def get_port():
        status, port = lpt.run_command('GetPort')
        response = make_response(port, status)
        response.mimetype = "text/plain"
        return response

    @app.put("/port")
    def set_port():
        port = request.get_data().decode("utf-8")
        status = lpt.run_command('SetPort',port)
        return port, status

    @app.get("/position")
    def get_position():
        status, pos = lpt.run_command('GetPosition')
        response = make_response(str(pos), status)
        response.mimetype = "text/plain"
        return response

    @app.put("/position")
    def set_position():
        pos = int(request.get_data().decode("utf-8"))
        status = lpt.run_command('SetPosition',pos)
        return str(pos), status

    return app

if __name__ == '__main__':
    app.run(threaded=False, processes=1)

# invoke from the command line as, for example:
#
# flask --app 'PYMEcs.Acquire.Hardware.NikonTiFlaskServerT:create_app(mode="production")' --debug run
