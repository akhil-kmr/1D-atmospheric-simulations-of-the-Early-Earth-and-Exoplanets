# stars.py

STARS = {
    "Sun": {
        "name": "Sun",
        "type": "G",
        "radius": 1.0,
        "temperature": 5778,
        "metallicity":0.0196,
        "distance": 1,     # pc
        "velocity": 0.0,     # km/s
        "planet name": "Earth",
        "planet radius": 6371,
        "planet diameter": 12742,
        "planet gravity": 9.81,
        "period": "365",
        "semi major axis": 1,
        "eccentricity": 0.0167,
        "stellar_frac": str(((12742e3/2)/(696340e3))**2)
    },
    
    "TRAPPIST-1": {
        "name": "TRAPPIST-1",
        "type": "M",
        "radius": 0.1192,
        "temperature": 2566,
        "metallicity":0.04,
        
        "planet name": "TRAPPIST-1 e",
        "planet radius": 6371,
        "planet diameter": 12742,
        "planet gravity": 9.81,
        "period": "365",
        "semi major axis": 1,
        "eccentricity": 0.0167,
        "stellar_frac": str(((12742e3/2)/(696340e3))**2)
    },
    
    "Proxima Centauri": {
        "name": "Proxima Centauri",
        "type": "M",
        "radius": 0.1542,
        "temperature": 2992,
        "metallicity":0.21,
        
        "planet name": "Proxima Centauri b",
        "planet radius": 6371,
        "planet diameter": 12742,
        "planet gravity": 9.81,
        "period": "365",
        "semi major axis": 1,
        "eccentricity": 0.0167,
        "stellar_frac": str(((12742e3/2)/(696340e3))**2)
    },
    
    "K2-18": {
        "name": "K2-18",
        "type": "M",
        "radius": 0.4445,
        "temperature": 3503,
        "metallicity":0.123,
        
        "planet name": "K2-18 b",
        "planet radius": 6371,
        "planet diameter": 12742,
        "planet gravity": 9.81,
        "period": "365",
        "semi major axis": 1,
        "eccentricity": 0.0167,
        "stellar_frac": str(((12742e3/2)/(696340e3))**2)
    },

    # Add more targets here
}

'''
elif (Star == 'K2V'):
    #Values for K2V star
    stellar_type = 'K'
    stellar_temp = '	4977'
    stellar_radius = '0.716'
    Metallicity = '−0.31'
    period = '30'
    eccentricity = '0'
    diameter = '12742' #planet diameter
    gravity = '9.81'   
    semi_major_axis = str(0.5) 
    stellar_frac = str(((12742e3/2)/(0.716*696340e3))**2)
        
elif (Star == 'TOI-270'):
    #Values for TOI-270 star
    stellar_type = 'M'
    stellar_temp = '3506'
    stellar_radius = '0.378'
    Metallicity = '-0.20'
    period = '11.379573'
    eccentricity = '0.0'
    diameter = str(int(12742*2.133)) #planet diameter
    gravity = '9.06'   
    total_T = round(total_T, 3)
    planet_distance = 22.477 #in parsec
    semi_major_axis = str(0.07210) # semi-major axis in AU
    stellar_frac = str(((12742e3*2.133/2)/(0.378*696340e3))**2)
    
elif (Star == 'TOI-776'):
    #Values for TOI-776 star
    stellar_type = 'M'
    stellar_temp = '3725'
    stellar_radius = '0.547'
    Metallicity = '-0.21'
    period = '15.665323'
    eccentricity = '0.0'
    diameter = str(int(12742*2.047)) #planet diameter
    gravity = '16.15'   
    total_T = round(total_T, 3)
    planet_distance = 27.29 #in parsec
    semi_major_axis = str(0.1001) # semi-major axis in AU
    stellar_frac = str(((12742e3*2.047/2)/(0.547*696340e3))**2)
    
elif (Star == 'K2-18'):
    #Values for K218 star
    stellar_type = 'M'
    stellar_temp = '3503'
    stellar_radius = '' #Benneke et al 2019
    Metallicity = '0.123'
    period = '32.94' #in days
    eccentricity = '0.0'
    diameter = str(int(12742*2.61)) #planet diameter
    gravity = '12.2'   
    total_T = round(total_T, 3)
    planet_distance = 38 #in parsec
    semi_major_axis = str(0.15910) # semi-major axis in AU
    stellar_frac = str(((12742e3*2.61/2)/(0.4445*696340e3))**2)
if (Star == 'Sun special'):
    #Values for Earth
    stellar_type = 'G' #stellar type
    stellar_temp = '5777' #stellar effective temperature
    stellar_radius = '1' #radius of star
    Metallicity = '0.0' #metallicity of star
    period = '365' #orbital period
    eccentricity = '0.01670' #orbital eccentricity
    diameter = str(int(12742*2.61)) #planet diameter
    gravity = '12.2' #gravity on the planet
    total_T = round(total_T, 3)
    planet_distance = 10 #in parsec
    semi_major_axis = str(semi_major) # semi-major axis in AU
    stellar_frac = str(((12742e3*2.61/2)/(696340e3))**2)
'''