from InitialPositions.randomUniform import create_randomUniform

def getInit(config):
    """Return function which returns initial positions

    Arguments:
        config {dictionary} -- Dictionary with parameters on how to do the initialization

    Returns:
        function -- A function which returns initial positions, made according to the config
    """
    if config["initial"]["name"] == "randomUniform":
        return create_randomUniform(config)
    else:
        raise NotImplementedError