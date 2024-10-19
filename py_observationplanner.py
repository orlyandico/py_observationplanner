import csv
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
import astropy.units as u
from astropy.coordinates import Angle
import numpy as np
from datetime import datetime, timedelta
from zoneinfo import ZoneInfo
from astroplan import Observer

def get_twilight_times(observer, start_time, end_time):
    # Convert to observer's timezone
    local_start = start_time.datetime.astimezone(observer.timezone)
    local_end = end_time.datetime.astimezone(observer.timezone)

    # If it's after midnight, we need to look at the previous evening for twilight start
    if local_start.hour < 12:
        evening_base = Time((local_start - timedelta(days=1)).replace(hour=12, minute=0, second=0, microsecond=0))
    else:
        evening_base = Time(local_start.replace(hour=12, minute=0, second=0, microsecond=0))

    night_start = observer.twilight_evening_astronomical(evening_base, which='next')
    night_end = observer.twilight_morning_astronomical(night_start, which='next')

    # Convert night_start and night_end to timezone-aware datetimes
    night_start_local = night_start.datetime.replace(tzinfo=observer.timezone)
    night_end_local = night_end.datetime.replace(tzinfo=observer.timezone) if night_end else None

    # If twilight end is after our end_time, it's not visible in our window
    if night_end_local and night_end_local > local_end:
        night_end_local = None

    # If twilight start is after our end_time, the entire night is not in our window
    if night_start_local > local_end:
        return None, None

    return night_start_local, night_end_local

def parse_ra_dec(ra_str, dec_str):
    ra = Angle(ra_str + ' hours')
    dec = Angle(dec_str + ' degrees')
    return SkyCoord(ra=ra, dec=dec, frame='icrs')

def read_csv(filename):
    objects = []
    with open(filename, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['RA'] and row['DEC']:
                coord = parse_ra_dec(row['RA'], row['DEC'])
                objects.append({
                    'ID': row['ID'],
                    'coord': coord,
                    'name': f"{row['ID']} ({row['Type']} in {row['Constellation']}) - {row['Size']}\n{row['Remarks']}"
                })
    return objects

def is_dark(time, observer):
    sun_altaz = get_sun(time).transform_to(AltAz(obstime=time, location=observer.location))
    return sun_altaz.alt < -18*u.deg  # Sun is 18 degrees below horizon (astronomical twilight)

def find_rise_set_times(coord, observer, times):
    altaz_frames = AltAz(obstime=[Time(t) for t in times], location=observer.location)
    obj_altazs = coord.transform_to(altaz_frames)

    above_horizon = obj_altazs.alt > 0*u.deg
    rises = np.where((~above_horizon[:-1]) & above_horizon[1:])[0]
    sets = np.where(above_horizon[:-1] & (~above_horizon[1:]))[0]

    rise_time = times[rises[0]] if len(rises) > 0 else None
    set_time = times[sets[0]] if len(sets) > 0 else None

    return rise_time, set_time

def find_window_times(coord, observer, alt, az, width, height, times):
    altaz_frames = AltAz(obstime=[Time(t) for t in times], location=observer.location)
    obj_altazs = coord.transform_to(altaz_frames)

    in_window = [is_in_fov(obj_altaz, alt, az, width, height) for obj_altaz in obj_altazs]
    entries = np.where((~np.array(in_window[:-1])) & np.array(in_window[1:]))[0]
    exits = np.where(np.array(in_window[:-1]) & (~np.array(in_window[1:])))[0]

    entry_time = times[entries[0]] if len(entries) > 0 else None
    exit_time = times[exits[0]] if len(exits) > 0 else None

    return entry_time, exit_time

def is_in_fov(obj_altaz, center_alt, center_az, width, height):
    alt_min = center_alt - height/2
    alt_max = center_alt + height/2
    az_min = (center_az - width/2) % 360
    az_max = (center_az + width/2) % 360

    obj_alt = obj_altaz.alt.deg
    obj_az = obj_altaz.az.deg

    in_alt_range = alt_min <= obj_alt <= alt_max

    if az_min < az_max:
        in_az_range = az_min <= obj_az <= az_max
    else:  # Field of view crosses 0/360 degree azimuth
        in_az_range = obj_az >= az_min or obj_az <= az_max

    return in_alt_range and in_az_range


def calculate_visibility(objects, observer, alt, az, width, height, start_time, end_time, timezone):
    results = []
    times = [start_time.datetime.replace(tzinfo=timezone) + timedelta(minutes=i*5) for i in range(int((end_time.datetime - start_time.datetime).total_seconds() / 300))]

    for obj in objects:
        rise_time, set_time = find_rise_set_times(obj['coord'], observer, times)
        entry_time, exit_time = find_window_times(obj['coord'], observer, alt, az, width, height, times)

        if rise_time or set_time or entry_time or exit_time:
            results.append({
                'name': obj['name'],
                'coord': obj['coord'],
                'rise_time': rise_time,
                'set_time': set_time,
                'entry_time': entry_time,
                'exit_time': exit_time,
                'has_dark_event': (rise_time and is_dark(Time(rise_time), observer)) or
                                  (set_time and is_dark(Time(set_time), observer)) or
                                  (entry_time and is_dark(Time(entry_time), observer)) or
                                  (exit_time and is_dark(Time(exit_time), observer))
            })

    return results


def get_altaz_at_time(coord, observer, time):
    altaz_frame = AltAz(obstime=Time(time), location=observer.location)
    altaz = coord.transform_to(altaz_frame)
    return altaz.alt.deg, altaz.az.deg

def main():
    csv_file = 'astronomical_objects.csv'  # Replace with your CSV file path
    lat = 51.5074  # London latitude
    lon = -0.1278  # London longitude
    alt = 60  # 45 degrees altitude
    az = 250  # 220 degrees azimuth
    width = 160  # 120 degrees width
    height = 60  # 60 degrees height

    # Set timezone for London
    timezone = ZoneInfo('Europe/London')

    # Create observer
    location = EarthLocation(lat=lat*u.deg, lon=lon*u.deg)
    observer = Observer(location=location, timezone=timezone)

    # Use current time as start time in the specified timezone
    current_time = datetime.now(timezone)
    start_time = Time(current_time)
    # Set end time to 24 hours from now
    end_time = Time(start_time.datetime + timedelta(days=1))

    # Get astronomical twilight times
    twilight_start, twilight_end = get_twilight_times(observer, start_time, end_time)

    objects = read_csv(csv_file)
    visible_objects = calculate_visibility(objects, observer, alt, az, width, height, start_time, end_time, timezone)

    # Sort objects based on the closest event (entry or exit) to current time
    def sort_key(obj):
        entry_diff = abs((obj['entry_time'] - current_time).total_seconds()) if obj['entry_time'] else float('inf')
        exit_diff = abs((obj['exit_time'] - current_time).total_seconds()) if obj['exit_time'] else float('inf')
        return min(entry_diff, exit_diff)

    sorted_objects = sorted(visible_objects, key=sort_key)

    print(f"Report for night of {start_time.datetime.strftime('%Y-%m-%d')}")
    print(f"Current time: {current_time.strftime('%Y-%m-%d %H:%M:%S %Z')}")
    print(f"Astronomical Twilight:")
    if twilight_start and twilight_end:
        print(f"  Start: {twilight_start.strftime('%Y-%m-%d %H:%M:%S %Z')}")
        print(f"  End: {twilight_end.strftime('%Y-%m-%d %H:%M:%S %Z')}")
    elif twilight_start:
        print(f"  Start: {twilight_start.strftime('%Y-%m-%d %H:%M:%S %Z')}")
        print("  End: After the end of the observation period")
    else:
        print("  No astronomical night during the observation period")
    print()

    print(f"Objects visible from {start_time.datetime.strftime('%Y-%m-%d %H:%M:%S %Z')} to {end_time.datetime.strftime('%Y-%m-%d %H:%M:%S %Z')}:")
    print(f"Window: {width}°x{height}° field centered at alt={alt}°, az={az}°")
    print("(* indicates event occurs during dark hours)")
    print("Sorted by closeness of entry or exit to current time:")
    print()

    for obj in sorted_objects:
        if obj['has_dark_event']:
            print(f"{obj['name']}:")
            if obj['rise_time']:
                dark_indicator = '*' if is_dark(Time(obj['rise_time']), observer) else ''
                print(f"  Rises: {obj['rise_time'].strftime('%Y-%m-%d %H:%M:%S %Z')}{dark_indicator}")
            if obj['set_time']:
                dark_indicator = '*' if is_dark(Time(obj['set_time']), observer) else ''
                print(f"  Sets: {obj['set_time'].strftime('%Y-%m-%d %H:%M:%S %Z')}{dark_indicator}")
            if obj['entry_time']:
                dark_indicator = '*' if is_dark(Time(obj['entry_time']), observer) else ''
                time_diff = obj['entry_time'] - current_time
                entry_alt, entry_az = get_altaz_at_time(obj['coord'], observer, obj['entry_time'])
                print(f"  Enters window: {obj['entry_time'].strftime('%Y-%m-%d %H:%M:%S %Z')}{dark_indicator}")
                print(f"    Alt/Az at entry: {entry_alt:.2f}°/{entry_az:.2f}°")
            if obj['exit_time']:
                dark_indicator = '*' if is_dark(Time(obj['exit_time']), observer) else ''
                time_diff = obj['exit_time'] - current_time
                print(f"  Exits window: {obj['exit_time'].strftime('%Y-%m-%d %H:%M:%S %Z')}{dark_indicator}")
            print()

if __name__ == "__main__":
    main()
