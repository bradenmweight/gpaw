# A. Pletzer Tue Mar 20 11:42:05 EST 2001
# G. Genellina 2009-09-10: Minor syntax changes,
#   compatibility with Python 2.4 and above
# downloaded from http://code.activestate.com/recipes/52292/
# and slightly changed

"""
Gauss-Legendre Integration
"""

from itertools import izip
import numpy as np

_nodes = (
    (0.,),
    (-0.5773502691896257,
      0.5773502691896257,),
    (-0.7745966692414834,
      0.,
      0.7745966692414834,),
    (-0.861136311594053,
      -0.3399810435848562,
      0.3399810435848562,
      0.861136311594053,),
    (-0.906179845938664,
      -0.5384693101056829,
      0.,
      0.5384693101056829,
      0.906179845938664,),
    (-0.932469514203152,
      -0.6612093864662646,
      -0.2386191860831968,
      0.2386191860831968,
      0.6612093864662646,
      0.932469514203152,),
    (-0.949107912342759,
      -0.7415311855993937,
      -0.4058451513773972,
      0.,
      0.4058451513773971,
      0.7415311855993945,
      0.949107912342759,),
    (-0.960289856497537,
      -0.7966664774136262,
      -0.5255324099163289,
      -0.1834346424956498,
      0.1834346424956498,
      0.5255324099163289,
      0.7966664774136262,
      0.960289856497537,),
    (-0.968160239507626,
      -0.836031107326637,
      -0.6133714327005903,
      -0.3242534234038088,
      0.,
      0.3242534234038088,
      0.6133714327005908,
      0.836031107326635,
      0.968160239507627,),
    (-0.973906528517172,
      -0.865063366688984,
      -0.6794095682990246,
      -0.433395394129247,
      -0.1488743389816312,
      0.1488743389816312,
      0.433395394129247,
      0.6794095682990246,
      0.865063366688984,
      0.973906528517172,),
    (-0.97822865814604,
      -0.88706259976812,
      -0.7301520055740422,
      -0.5190961292068116,
      -0.2695431559523449,
      0.,
      0.2695431559523449,
      0.5190961292068117,
      0.73015200557405,
      0.887062599768093,
      0.978228658146058,),
    (-0.981560634246732,
      -0.904117256370452,
      -0.7699026741943177,
      -0.5873179542866143,
      -0.3678314989981804,
      -0.1252334085114688,
      0.1252334085114688,
      0.3678314989981804,
      0.5873179542866143,
      0.7699026741943177,
      0.904117256370452,
      0.981560634246732,),
    )

_weights = (
    (2.,),
    (1.,
     1.,),
    (0.5555555555555553,
     0.888888888888889,
     0.5555555555555553,),
    (0.3478548451374539,
     0.6521451548625462,
     0.6521451548625462,
     0.3478548451374539,),
    (0.2369268850561887,
     0.4786286704993665,
     0.5688888888888889,
     0.4786286704993665,
     0.2369268850561887,),
    (0.1713244923791709,
     0.3607615730481379,
     0.4679139345726913,
     0.4679139345726913,
     0.3607615730481379,
     0.1713244923791709,),
    (0.129484966168868,
     0.2797053914892783,
     0.3818300505051186,
     0.4179591836734694,
     0.3818300505051188,
     0.279705391489276,
     0.1294849661688697,),
    (0.1012285362903738,
     0.2223810344533786,
     0.3137066458778874,
     0.3626837833783619,
     0.3626837833783619,
     0.3137066458778874,
     0.2223810344533786,
     0.1012285362903738,),
    (0.0812743883615759,
     0.1806481606948543,
     0.2606106964029356,
     0.3123470770400029,
     0.3302393550012597,
     0.3123470770400025,
     0.2606106964029353,
     0.1806481606948577,
     0.0812743883615721,),
    (0.06667134430868681,
     0.149451349150573,
     0.2190863625159832,
     0.2692667193099968,
     0.2955242247147529,
     0.2955242247147529,
     0.2692667193099968,
     0.2190863625159832,
     0.149451349150573,
     0.06667134430868681,),
    (0.05566856711621584,
     0.1255803694648743,
     0.1862902109277404,
     0.2331937645919927,
     0.2628045445102466,
     0.2729250867779006,
     0.2628045445102466,
     0.2331937645919933,
     0.1862902109277339,
     0.1255803694649132,
     0.05566856711616958,),
    (0.04717533638647547,
     0.1069393259953637,
     0.1600783285433586,
     0.2031674267230672,
     0.2334925365383534,
     0.2491470458134027,
     0.2491470458134027,
     0.2334925365383534,
     0.2031674267230672,
     0.1600783285433586,
     0.1069393259953637,
     0.04717533638647547,),
    )

_nodesLog = (
    (0.3333333333333333,),
    (0.1120088061669761,
     0.6022769081187381,),
    (0.06389079308732544,
     0.3689970637156184,
     0.766880303938942,),
    (0.04144848019938324,
     0.2452749143206022,
     0.5561654535602751,
     0.848982394532986,),
    (0.02913447215197205,
     0.1739772133208974,
     0.4117025202849029,
     0.6773141745828183,
     0.89477136103101,),
    (0.02163400584411693,
     0.1295833911549506,
     0.3140204499147661,
     0.5386572173517997,
     0.7569153373774084,
     0.922668851372116,),
    (0.0167193554082585,
     0.100185677915675,
     0.2462942462079286,
     0.4334634932570557,
     0.6323509880476823,
     0.81111862674023,
     0.940848166743287,),
    (0.01332024416089244,
     0.07975042901389491,
     0.1978710293261864,
     0.354153994351925,
     0.5294585752348643,
     0.7018145299391673,
     0.849379320441094,
     0.953326450056343,),
    (0.01086933608417545,
     0.06498366633800794,
     0.1622293980238825,
     0.2937499039716641,
     0.4466318819056009,
     0.6054816627755208,
     0.7541101371585467,
     0.877265828834263,
     0.96225055941096,),
    (0.00904263096219963,
     0.05397126622250072,
     0.1353118246392511,
     0.2470524162871565,
     0.3802125396092744,
     0.5237923179723384,
     0.6657752055148032,
     0.7941904160147613,
     0.898161091216429,
     0.9688479887196,),
    (0.007643941174637681,
     0.04554182825657903,
     0.1145222974551244,
     0.2103785812270227,
     0.3266955532217897,
     0.4554532469286375,
     0.5876483563573721,
     0.7139638500230458,
     0.825453217777127,
     0.914193921640008,
     0.973860256264123,),
    (0.006548722279080035,
     0.03894680956045022,
     0.0981502631060046,
     0.1811385815906331,
     0.2832200676673157,
     0.398434435164983,
     0.5199526267791299,
     0.6405109167754819,
     0.7528650118926111,
     0.850240024421055,
     0.926749682988251,
     0.977756129778486,),
    )

_weightsLog = (
    (-1.,),
    (-0.7185393190303845,
      -0.2814606809696154,),
    (-0.5134045522323634,
      -0.3919800412014877,
      -0.0946154065661483,),
    (-0.3834640681451353,
      -0.3868753177747627,
      -0.1904351269501432,
      -0.03922548712995894,),
    (-0.2978934717828955,
      -0.3497762265132236,
      -0.234488290044052,
      -0.0989304595166356,
      -0.01891155214319462,),
    (-0.2387636625785478,
      -0.3082865732739458,
      -0.2453174265632108,
      -0.1420087565664786,
      -0.05545462232488041,
      -0.01016895869293513,),
    (-0.1961693894252476,
      -0.2703026442472726,
      -0.239681873007687,
      -0.1657757748104267,
      -0.0889432271377365,
      -0.03319430435645653,
      -0.005932787015162054,),
    (-0.164416604728002,
      -0.2375256100233057,
      -0.2268419844319134,
      -0.1757540790060772,
      -0.1129240302467932,
      -0.05787221071771947,
      -0.02097907374214317,
      -0.003686407104036044,),
    (-0.1400684387481339,
      -0.2097722052010308,
      -0.211427149896601,
      -0.1771562339380667,
      -0.1277992280331758,
      -0.07847890261203835,
      -0.0390225049841783,
      -0.01386729555074604,
      -0.002408041036090773,),
    (-0.12095513195457,
      -0.1863635425640733,
      -0.1956608732777627,
      -0.1735771421828997,
      -0.135695672995467,
      -0.0936467585378491,
      -0.05578772735275126,
      -0.02715981089692378,
      -0.00951518260454442,
      -0.001638157633217673,),
    (-0.1056522560990997,
      -0.1665716806006314,
      -0.1805632182877528,
      -0.1672787367737502,
      -0.1386970574017174,
      -0.1038334333650771,
      -0.06953669788988512,
      -0.04054160079499477,
      -0.01943540249522013,
      -0.006737429326043388,
      -0.001152486965101561,),
    (-0.09319269144393,
      -0.1497518275763289,
      -0.166557454364573,
      -0.1596335594369941,
      -0.1384248318647479,
      -0.1100165706360573,
      -0.07996182177673273,
      -0.0524069547809709,
      -0.03007108900074863,
      -0.01424924540252916,
      -0.004899924710875609,
      -0.000834029009809656,),
    )

_NGMAX = len(_nodes)
_NGMIN = 1

class GaussLegendre:
    def __init__(self, xmin, xmax, ng=10):
        try:
            ns = _nodes[ng-1]
            ws = _weights[ng-1]
        except:
            raise RuntimeError, ('Gauss-Legendre only possible for n=' +
                                 str(_NGMIN) + '-' + str(_NGMAX))
        dx = xmax - xmin
        self.x = np.array([(dx * y + xmin + xmax)/2. for y in ns])
        self.w = np.array([(0.5 * dx * w) for w in ws])

    def get_x(self):
        """Return abscissas."""
        return self.x

    def get_w(self):
        """Return weights."""
        return self.w

def gauss(xmin, xmax, funct, ng=10):
    """
    Gauss quadature (weight function = 1.0):
    xmin, xmax: boundaries of integration domain
    funct: integrand function
    ng: Gauss integration order
    """
    ng = max(min(ng, _NGMAX), _NGMIN)
    ns = _nodes[ng-1]
    ws = _weights[ng-1]
    dx = xmax - xmin
    xs = [(dx*y + xmin + xmax)/2. for y in ns]
    return 0.5*dx*sum(funct(x)*w for x,w in izip(xs,ws))

def gaussLog(xmin, xmax, funct, ng=10):
    """
    Gauss quadature with Log singularity at x=xmin:
    xmin, xmax: boundaries of integration domain
    funct: integrand function
    ng: Gauss integration order
    """
    ng = max(min(ng, _NGMAX), _NGMIN)
    ns = _nodesLog[ng-1]
    ws = _weightsLog[ng-1];
    dx = xmax - xmin
    xs = [(dx*y + xmin) for y in ns]
    return dx*sum(funct(x)*w for x,w in izip(xs,ws))

####

if __name__ == '__main__':

    from math import *

    def f2(x):
        return x**2

    def f3(x):
        return x**4

    def f4(x):
        return cos(2.*pi*(x-0.128726465))

    def f5(x):
        return 2.*cos(2.*pi*(x-0.128726465))**2

    print('-'*80)
    print('Gauss (weight function = 1)')
    print('-'*80)

    # simple tests
    print('gauss(0., 1., f3, 1)=', gauss(0., 1., f3, 1))
    print('gauss(0., 1., f3, 2)=', gauss(0., 1., f3, 2))
    print('gauss(0., 1., f3, 3)=', gauss(0., 1., f3, 3))
    print('gauss(0., 1., f3   )=', gauss(0., 1., f3   ))

    # convergence test
    ng = range(_NGMIN, _NGMAX+1)

    print("""\n
    Integrate[Cos[2.*Pi*(x-0.128726465)], {x, 0, 1}]
    \n""")

    error = [gauss(0., 10.0, f4, n) for n in ng]
    print('    n = ', '%8d'*len(ng) % tuple(ng))
    print('error = ', '%8.0e'*len(error) % tuple(error))

    print("""\n
    Integrate[2.*Cos[2.*Pi*(x-0.128726465)]^2, {x, 0, 1}]
    \n""")

    error = [gauss(0., 1.0, f5, n)-1.0 for n in ng]
    print('    n = ', '%8d'*len(ng) % tuple(ng))
    print('error = ', '%8.0e'*len(error) % tuple(error))


    print('-'*80)
    print('Gauss with Log singularity at left boundary')
    print('-'*80)

    a, b = 0., 1.

    print("""\n
    Integrate[Log[x]*x^2, {x, 0, 1}]
    \n""")

    exact = -1./9.
    error = [gaussLog(a, b, f2, n) - exact for n in ng]
    print('    n = ', '%8d'*len(ng) % tuple(ng))
    print('error = ', '%8.0e'*len(error) % tuple(error))

    print("""\n
    Integrate[Log[x]*2.*Cos[2.*Pi*(x-0.128726465)]^2, {x, 0, 1}]
    \n""")

    exact = -1.242002481967963
    error = [gaussLog(a, b, f5, n) - exact for n in ng]
    print('    n = ', '%8d'*len(ng) % tuple(ng))
    print('error = ', '%8.0e'*len(error) % tuple(error))
