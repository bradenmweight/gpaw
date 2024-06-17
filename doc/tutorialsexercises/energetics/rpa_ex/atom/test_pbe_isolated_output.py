import numpy as np


def test():
    isolated_results, isolated_benchmark = [], [
        [6.0, -0.664402266578, 9.88627820237],
        [7.0, -0.778484948334, 9.79998115788],
        [8.0, -0.82500272946, 9.76744817185],
        [9.0, -0.841856681349, 9.75715732758],
        [10.0, -0.848092042293, 9.75399390142],
        [11.0, -0.850367362642, 9.75296805021],
        [12.0, -0.85109735188, 9.75265131464]
    ]
    with open("si.atom.pbe_and_exx_energies.txt", "r") as file:
        lines = file.readlines()
        for i, line in enumerate(lines):
            if i == 0:
                continue
            x = [float(x) for x in line.strip().split()]
            isolated_results.append(x)

    for result, benchmark in zip(isolated_results, isolated_benchmark):
        assert len(result) == 3
        print(result, benchmark)
        assert np.allclose(result, benchmark, rtol=1.e-5, atol=1.e-8)


if __name__ == "__main__":
    test()
