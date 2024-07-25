### path setup ####################################################################################

import datetime

###################################################################################################
### Project Information ###########################################################################
###################################################################################################

project = 'dynamic_characterization'
author = 'Timo Diepers'
copyright = datetime.date.today().strftime("%Y") + ' brightway developers'
version: str = 'latest' # required by the version switcher

###################################################################################################
### Project Configuration #########################################################################
###################################################################################################

needs_sphinx = '7.3.0'

extensions = [
    # core extensions
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',    
    'sphinx.ext.intersphinx',
    'sphinx.ext.extlinks',
    'sphinx.ext.inheritance_diagram',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    # iPython extensions
    'IPython.sphinxext.ipython_directive',
    'IPython.sphinxext.ipython_console_highlighting',
    # Markdown support
    # 'myst_parser', # do not enable separately if using myst_nb, compare: https://github.com/executablebooks/MyST-NB/issues/421#issuecomment-1164427544
    # Jupyter Notebook support
    'myst_nb',
    # API documentation support
    'autoapi',
    # responsive web component support
    'sphinx_design',
    # custom 404 page
    'notfound.extension',
    # custom favicons
    'sphinx_favicon',
    # copy button on code blocks
    "sphinx_copybutton",
]

autoapi_dirs = ['../dynamic_characterization']
autoapi_type = 'python'
autoapi_ignore = [
    '*/data/*',
    '*tests/*',
    '*tests.py',
    '*validation.py',
    '*version.py',
    '*.rst',
    '*.yml',
    '*.md',
    '*.json',
    '*.data'
]

autoapi_options = [
    'members',
    'undoc-members',
    'private-members',
    'show-inheritance',
    'show-module-summary',
    #'special-members',
    #'imported-members',
    'show-inheritance-diagram'
]

autoapi_python_class_content = 'both'
autoapi_member_order = 'bysource'
autoapi_root = 'content/api'
autoapi_keep_files = False


autosummary_generate = True

master_doc = "index"

root_doc = 'index'
html_static_path = ['_static']
templates_path = ['_templates']
exclude_patterns = ['_build']
html_theme = "pydata_sphinx_theme"

suppress_warnings = [
    "myst.header" # suppress warnings of the kind "WARNING: Non-consecutive header level increase; H1 to H3"
]


####################################################################################################
### Theme html Configuration #######################################################################
####################################################################################################

html_show_sphinx = False
html_show_copyright = True

html_css_files = [
    "custom.css",
    "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.1.1/css/all.min.css" # for https://fontawesome.com/ icons
]

html_sidebars = {
    "**": [
        "sidebar-nav-bs.html",
    ],
    "content/index": [],
    "content/installation": [],
    "content/theory": [],
    "content/contributing": [],
    "content/codeofconduct": [],
    "content/license": [],
    "content/chagnelog": [],
}

html_theme_options = {
    # page elements
    "announcement": "⚠️ This package is under active development and some functionalities may change in the future.",
    "navbar_start": ["navbar-logo"],
    "navbar_end": ["theme-switcher", "navbar-icon-links.html"],
    "navbar_align": "left",
    # "navbar_persistent": ["theme-switcher"], # this is where the search button is usually placed
    "footer_start": ["copyright"],
    "footer_end": ["footer"],
    "secondary_sidebar_items": ["page-toc", "edit-this-page", "sourcelink", "support"],
    "header_links_before_dropdown": 7,
    # page elements content
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/TimoDiepers/dynamic_characterization",
            "icon": "fab fa-brands fa-github",
        },
        # {
        #     "name": "Conda",
        #     "url": "https://anaconda.org/diepers/dynamic_characterization",
        #     "icon": "fa-brands fa-python",
        #     "type": "fontawesome",
        # },
    ],
    # various settings
    "collapse_navigation": True,
    "show_prev_next": False,
    "use_edit_page_button": True,
    "navigation_with_keys": True,
    "logo": {
        "text": "dynamic_characterization",
        "image_light": "https://raw.githubusercontent.com/brightway-lca/brightway-documentation/main/source/_static/logo/BW_all_black_transparent_landscape.svg",
        "image_dark": "https://raw.githubusercontent.com/brightway-lca/brightway-documentation/main/source/_static/logo/BW_all_white_transparent_landscape.svg"
    },
}

# required by html_theme_options: "use_edit_page_button"
html_context = {
    "github_user": "TimoDiepers",
    "github_repo": "dynamic_characterization",
    "github_version": "main",
    "doc_path": "docs",
}

# notfound Configuration ################################################
# https://sphinx-notfound-page.readthedocs.io

notfound_context = {
    'title': 'Page Not Found',
    'body': '''                                                                                                                                           
        <h1>🍂 Page Not Found (404)</h1>
        <p>
        Oops! It looks like you've stumbled upon a page that's been recycled into the digital abyss.
        But don't worry, we're all about sustainability here.
        Why not take a moment to reduce, reuse, and recycle your clicks by heading back to the main page?
        And remember, every little bit counts in the grand scheme of things.
        </p>
    ''',
}

####################################################################################################
### Extension Configuration ########################################################################
####################################################################################################

# myst_parser Configuration ############################################
# https://myst-parser.readthedocs.io/en/latest/configuration.html

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'myst-nb',
    '.ipynb': 'myst-nb'
}


myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_image",
]

# myst-nb configuration ################################################
# https://myst-nb.readthedocs.io/en/latest/configuration.html

nb_execution_mode = 'off'

# sphinx-favicon configuration #########################################
# https://github.com/tcmetzger/sphinx-favicon

favicons = [
    {
        "rel": "icon", 
        "href": "favicon.svg", 
        "type": "image/svg+xml"
    },    
    {
        "rel": "icon", 
        "sizes": "144x144",
        "href": "favicon-144.png", 
        "type": "image/png"
    },
    {
        "rel": "mask-icon",
        "href": "favicon_mask-icon.svg",
        "color": "#222832"
    },
    {
        "rel": "apple-touch-icon",
        "sizes": "500x500",
        "href": "favicon-500.png"
    },
]
