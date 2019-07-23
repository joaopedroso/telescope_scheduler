"""
hidden: Indicate positions unavailable for observations

Function `hidden(sky)` is called before each observation to indicate positions not
currently available, typically due to being concealed by clouds
"""


def hidden(sky):
    """
    Return indices in `sky` list that are currently hidden by clouds or other objects

    .. note::
        `sky` list positions hidden by the moon are removed during preprocessing

    :param sky: list of `SkyCoord` positions to be observed
    :return: indices of `sky` that cannot be observed now
    :rtype: sequence of int
    """
    # return set(range(100,110))   # usage example
    return []
