# scripts

Loose collection of scripts for astronomy research.

1. besancon_model/besancon_gaia_dr2.py

    The output of python besancon_model/besancon_gaia_dr2.py --help is shown below:
    
    This Script downloads a besancon model including include_kinematics using the
astroquery interface. Gaia DR2 parallax and proper motion sky-averaged
uncertainties are derived using pygaia. Random errors corresponding to Gaia
DR2 expectations can be added to the parallaxes and distances of the model
output.

    positional arguments:
      email                 Your valid email address.
    
    optional arguments:
      -h, --help            show this help message and exit  
      --ra_deg RA_DEG       RA in degrees  
      --dec_deg DEC_DEG     Dec in degrees  
      --field_size_squaredegree FIELD_SIZE_SQUAREDEGREE  
                            Field size in squaredegrees.  
      --data_dir DATA_DIR   Path to directory for table saving  
      --overwrite OVERWRITE Overwrite download from besancon server.  
      --add_random_gaia_errors ADD_RANDOM_GAIA_ERRORS  
                            Add random errors to parallax and PM assuming Gaia DR2  
                            performances.  
      --random_seed RANDOM_SEED  
                            numpy random seed for error simulation. Allows for  
                            repeatability.  
