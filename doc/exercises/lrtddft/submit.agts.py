from q2.job import Job


def workflow():
    return [
        Job('Na2TDDFT.py@2x1m'),
        Job('part2.py', deps=['Na2TDDFT.py']),
        Job('ground_state.py@8x15s'),
        Job('spectrum.py', deps=['ground_state.py'])]
