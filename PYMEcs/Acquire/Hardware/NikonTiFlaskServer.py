from flask import Flask, request, jsonify, make_response
from http import HTTPStatus

def create_app(mode='simulated'):
    # create and configure the app
    app = Flask(__name__)
    
    if mode == 'simulated':
        from PYMEcs.Acquire.Hardware import NikonTiSim
        lp = NikonTiSim.LightPath()
    else:
        from PYME.Acquire.Hardware import NikonTi
        lp = NikonTi.LightPath()
        # lp.SetPosition(2)

    @app.get("/names")
    def get_names():
        return jsonify(lp.names)

    @app.get("/port")
    def get_port():
        response = make_response(lp.GetPort(), HTTPStatus.OK)
        response.mimetype = "text/plain"
        return response

    @app.put("/port")
    def set_port():
        lp.SetPort(request.get_data().decode("utf-8"))
        return lp.GetPort(), HTTPStatus.OK

    @app.get("/position")
    def get_position():
        response = make_response(str(lp.GetPosition()), HTTPStatus.OK)
        response.mimetype = "text/plain"
        return response

    @app.put("/position")
    def set_position():
        lp.SetPosition(int(request.get_data().decode("utf-8")))
        return str(lp.GetPosition()), HTTPStatus.OK

    return app

if __name__ == '__main__':
    app.run(threaded=False, processes=1)
