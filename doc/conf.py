import datetime
import sys

import sphinx_rtd_theme
from gpaw import __version__
from gpaw.doctools.aamath import autodoc_process_docstring
try:
    import sphinxcontrib.spelling
except ImportError:
    sphinxcontrib = None

assert sys.version_info >= (3, 8)

sys.path.append('.')

extensions = ['images',
              'ext',
              'sphinx.ext.autodoc',
              'sphinx.ext.doctest',
              'sphinx.ext.extlinks',
              'sphinx.ext.viewcode',
              'sphinx.ext.napoleon',
              'sphinx.ext.mathjax',
              'sphinx.ext.intersphinx']

if sphinxcontrib:
    extensions.append('sphinxcontrib.spelling')
extlinks = {
    'doi': ('https://doi.org/%s', 'doi: %s'),
    'arxiv': ('https://arxiv.org/abs/%s', 'arXiv: %s'),
    'mr': ('https://gitlab.com/gpaw/gpaw/-/merge_requests/%s', 'MR: !%s'),
    'issue': ('https://gitlab.com/gpaw/gpaw/-/issues/%s', 'issue: #%s'),
    'xkcd': ('https://xkcd.com/%s', 'XKCD: %s')}
spelling_word_list_filename = 'words.txt'
spelling_show_suggestions = True
templates_path = ['templates']
source_suffix = '.rst'
master_doc = 'index'
project = 'GPAW'
copyright = f'{datetime.date.today().year}, GPAW developers'
release = __version__
exclude_patterns = ['build']
default_role = 'math'
pygments_style = 'sphinx'
autoclass_content = 'both'
modindex_common_prefix = ['gpaw.']
intersphinx_mapping = {
    'python': ('https://docs.python.org/3.10', None),
    'ase': ('https://wiki.fysik.dtu.dk/ase', None),
    'numpy': ('https://numpy.org/doc/stable', None),
    'cupy': ('https://docs.cupy.dev/en/stable', None),
    'scipy': ('https://docs.scipy.org/doc/scipy', None),
    'pytest': ('https://docs.pytest.org/en/stable', None),
    'mayavi': ('https://docs.enthought.com/mayavi/mayavi', None)}
nitpick_ignore = [('py:class', 'gpaw.calculator.GPAW'),
                  ('py:class', 'gpaw.spinorbit.BZWaveFunctions'),
                  ('py:class', 'GPAW'),
                  ('py:class', 'Atoms'),
                  ('py:class', 'np.ndarray'),
                  ('py:class', 'ase.spectrum.dosdata.GridDOSData'),
                  ('py:class', 'ase.atoms.Atoms'),
                  ('py:class', 'gpaw.point_groups.group.PointGroup'),
                  ('py:class', 'UGArray'),
                  ('py:class', 'gpaw.core.arrays.DomainType'),
                  ('py:class', 'DomainType'),
                  ('py:class', 'Path'),
                  ('py:class', 'Vector'),
                  ('py:class', 'ArrayLike1D'),
                  ('py:class', 'ArrayLike2D'),
                  ('py:class', 'Array1D'),
                  ('py:class', 'Array2D'),
                  ('py:class', 'Array3D'),
                  ('py:class', 'MPIComm'),
                  ('py:class', 'NLOData'),
                  ('py:class', 'numpy._typing._dtype_like._SupportsDType'),
                  ('py:class', 'numpy._typing._dtype_like._DTypeDict'),
                  ('py:class', 'DTypeLike'),
                  ('py:class', 'ModuleType'),
                  ('py:class', 'IO'),
                  ('py:class', 'gpaw.gpu.cpupy.ndarray')]
html_theme = 'sphinx_rtd_theme'
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
html_style = 'gpaw.css'
html_title = 'GPAW'
html_favicon = 'static/gpaw_favicon.ico'
html_static_path = ['static']
html_last_updated_fmt = '%a, %d %b %Y %H:%M:%S'

mathjax3_config = {
    'tex': {
        'macros': {
            'br': '{\\mathbf r}',
            'bk': '{\\mathbf k}',
            'bG': '{\\mathbf G}'}}}

autodoc_typehints = 'description'
autodoc_typehints_description_target = 'documented'


def setup(app):
    app.connect('autodoc-process-docstring',
                lambda app, what, name, obj, options, lines:
                    autodoc_process_docstring(lines))
