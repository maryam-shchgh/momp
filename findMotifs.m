function minLocations = findMotifs(mp)

    minval = min(mp);
    minLocations = find(mp == minval);