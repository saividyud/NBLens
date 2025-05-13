
G = 6.6743e-20 # kg^-1 km^3 s^-2
AU2KM = 1.49597870700e8 # km

sun = {
    'name': 'Sun',
    'mass': 1.989e30, # kg
    'mu': 1.3271244e11, # km^3 / s^2
    'radius': 6.957e5, # km
    'colors': 'autumn'
}

earth = {
    'name': 'Earth',

    'mass': 5.972e24, # kg
    'mu': 3.98600436e5, # km^3 / s^2
    'radius': 6378.1, # km
    'a': 1.0000261, # au
    'e': 0.01671123,
    'i': -0.00078890, # degrees
    'P': 365.25636, # days

    'J2': 1.082635854e-3,

    'colors': 'winter_r'
}

moon = {
    'name': 'Moon',
    'mass': 0.07346e24, # kg
    'mu': 0.00490e6, # km^3 / s^2
    'radius': 1738.1, # km
    'a': 0.3844e6, # km
    'e': 0.0549,
    'i': 5.145, # degrees

    'colors': 'gray'
}

mercury = {
    'name': 'Mercury',
    'mass': 0.330103e24, # kg
    'mu': 0.022032e6, # km^3 / s^2
    'radius': 2439.4, # km
    'a': 0.38709927, # au
    'e': 0.20563593,
    'i': 7.00497902, # degrees
    'P': 74.8202572, # days
    'colors': 'gray'
}

venus = {
    'name': 'Venus',
    'mass': 4.86731e24, # kg
    'mu': 0.32486e6, # km^3 / s^2
    'radius': 6051.8, # km
    'a': 0.72333566, # au
    'e': 0.00677672,
    'i': 3.39467605, # degrees
    'P': 224.700799, # days
    'colors': 'Wistia'
}

mars = {
    'name': 'Mercury',
    'mass': 0.641691e24, # kg
    'mu': 42828.37362, # km^3 / s^2
    'radius': 3389.50, # km
    'a': 1.52371034, # au
    'e': 0.09339410,
    'i': 1.84969142, # degrees
    'P': 686.979586, # days
    'colors': 'Reds'
}

jupiter = {
    'name': 'Jupiter',
    'mass': 1898.12e24, # kg
    'mu': 126686531.9, # km^3 / s^2
    'radius': 69911, # km
    'a': 5.20288700, # au
    'e': 0.04838624,
    'i': 1.30439695, # degrees
    'P': 4332.820129, # days
    'colors': 'YlOrBr'
}

saturn = {
    'name': 'Saturn',
    'mass': 568.317e24, # kg
    'mu': 37931206.23, # km^3 / s^2
    'radius': 58232, # km
    'a': 9.53667594, # au
    'e': 0.05386179,
    'i': 2.48599187, # degrees
    'P': 10755.69864, # days
    'colors': 'Oranges'
}

uranus = {
    'name': 'Uranus',
    'mass': 86.8099e24, # kg
    'mu': 5793951.3, # km^3 / s^2
    'radius': 25362, # km
    'a': 19.18916464, # au
    'e': 0.04725744,
    'i': 0.77263783, # degrees
    'P': 30687.1530, # days
    'colors': 'Blues'
}

neptune = {
    'name': 'Neptune',
    'mass': 102.4092e24, # kg
    'mu': 6835099.97, # km^3 / s^2
    'radius': 24622, # km
    'a': 30.06992276, # au
    'e': 0.00859048,
    'i': 1.77004347, # degrees
    'P': 60190.0296, # days
    'colors': 'Blues'
}

pluto = {
    'name': 'Pluto',
    'mass': 13029e18, # kg
    'mu': 869.61, # km^3 / s^2
    'radius': 1188.3, # km
    'a': 39.4450697, # au
    'e': 0.25024871,
    'i': 17.0890009, # degrees
    'P': 90487.2769, # days
    'colors': 'pink'
}

ceres = {
    'name': 'Ceres',
    'mass': 938.416e18, # kg
    'mu': 938.416e18 * G, # km^3 / s^2
    'radius': 469.7, # km
    'a': 2.76661904, # au
    'e': 0.07863576,
    'i': 10.5867951, # degrees
    'P': 1680.824888, # days
    'colors': 'gray'
}

eris = {
    'name': 'Eris',
    'mass': 16600e18, # kg
    'mu': 16600e18 * G, # km^3 / s^2
    'radius': 1200, # km
    'a': 68.1175851, # au
    'e': 0.43196754,
    'i': 43.793657, # degrees
    'P': 205346.493, # days
    'colors': 'gist_yarg'
}

haumea = {
    'name': 'Haumea',
    'mass': 3100e18, # kg
    'mu': 3100e18 * G, # km^3 / s^2
    'radius': 714, # km
    'a': 42.9414319, # au
    'e': 0.19974147,
    'i': 28.2114985, # degrees
    'P': 102781.0880, # days
    'colors': 'binary'
}

makemake = {
    'name': 'Makemake',
    'mass': 4006, # kg
    'mu': 4006 * G, # km^3 / s^2
    'radius': 715, # km
    'a': 45.2593495, # au
    'e': 0.16607929,
    'i': 29.0101489, # degrees
    'P': 111214.4009, # days
    'colors': 'copper'
}

