def hidden(sky):
    """
    Return indices in `sky` list that are currently hidden by clouds or other objects

    .. note::
        `sky` list has positions hidden by the moon removed during preprocessing

    :param sky: list of `SkyCoord` positions to be observed
    :return: indices of `sky` that cannot be observed now
    :rtype: sequence of int
    """
    # return set(range(100,110))   # just as an example
    return []