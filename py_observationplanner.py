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
    local_start = start_time.datetime.astimezone(observer.timezone)
    local_end = end_time.datetime.astimezone(observer.timezone)

    evening_base = Time((local_start - timedelta(days=1 if local_start.hour < 12 else 0)).replace(hour=12, minute=0, second=0, microsecond=0))

    night_start = observer.twilight_evening_astronomical(evening_base, which='next')
    night_end = observer.twilight_morning_astronomical(night_start, which='next')

    night_start_local = night_start.datetime.replace(tzinfo=observer.timezone)
    night_end_local = night_end.datetime.replace(tzinfo=observer.timezone) if night_end else None

    if night_end_local and night_end_local > local_end:
        night_end_local = None

    return (None, None) if night_start_local > local_end else (night_start_local, night_end_local)

def parse_ra_dec(ra_str, dec_str):
    return SkyCoord(ra=Angle(ra_str + ' hours'), dec=Angle(dec_str + ' degrees'), frame='icrs')

def read_csv(filename):
    with open(filename, 'r') as f:
        reader = csv.DictReader(f)
        return [
            {
                'ID': row['ID'],
                'coord': parse_ra_dec(row['RA'], row['DEC']) if row['RA'] and row['DEC'] else None,
                'name': f"{row['ID']} ({row['Type']} in {row['Constellation']}) - {row['Size']}\n{row['Remarks']}"
            }
            for row in reader if row['RA'] and row['DEC']
        ]

def is_dark(time, observer):
    return get_sun(time).transform_to(AltAz(obstime=time, location=observer.location)).alt < -18*u.deg

def find_event_times(coord, observer, times, event_condition):
    altaz_frames = AltAz(obstime=[Time(t) for t in times], location=observer.location)
    obj_altazs = coord.transform_to(altaz_frames)

    events = np.where(event_condition(obj_altazs))[0]
    return times[events[0]] if len(events) > 0 else None

def is_in_fov(altazs, center_alt, center_az, width, height):
    alt_min, alt_max = center_alt - height/2, center_alt + height/2
    az_min, az_max = (center_az - width/2) % 360, (center_az + width/2) % 360

    obj_alts, obj_azs = altazs.alt.deg, altazs.az.deg

    in_alt_range = (alt_min <= obj_alts) & (obj_alts <= alt_max)

    if az_min < az_max:
        in_az_range = (az_min <= obj_azs) & (obj_azs <= az_max)
    else:  # Field of view crosses 0/360 degree azimuth
        in_az_range = (obj_azs >= az_min) | (obj_azs <= az_max)

    return in_alt_range & in_az_range

def find_event_times(coord, observer, times, event_condition):
    altaz_frames = AltAz(obstime=[Time(t) for t in times], location=observer.location)
    obj_altazs = coord.transform_to(altaz_frames)

    events = event_condition(obj_altazs)
    event_indices = np.where(events)[0]

    return [times[i] for i in event_indices] if len(event_indices) > 0 else []

def calculate_visibility(objects, observer, alt, az, width, height, start_time, end_time, timezone):
    times = [start_time.datetime.replace(tzinfo=timezone) + timedelta(minutes=i*5) for i in range(int((end_time.datetime - start_time.datetime).total_seconds() / 300))]

    results = []
    for obj in objects:
        altaz_frames = AltAz(obstime=[Time(t) for t in times], location=observer.location)
        obj_altazs = obj['coord'].transform_to(altaz_frames)

        above_horizon = obj_altazs.alt > 0*u.deg
        in_fov = is_in_fov(obj_altazs, alt, az, width, height)

        rise_times = find_event_times(obj['coord'], observer, times, lambda altazs: (~above_horizon[:-1]) & above_horizon[1:])
        set_times = find_event_times(obj['coord'], observer, times, lambda altazs: above_horizon[:-1] & (~above_horizon[1:]))
        entry_times = find_event_times(obj['coord'], observer, times, lambda altazs: (~in_fov[:-1]) & in_fov[1:])
        exit_times = find_event_times(obj['coord'], observer, times, lambda altazs: in_fov[:-1] & (~in_fov[1:]))

        if any([rise_times, set_times, entry_times, exit_times]):
            results.append({
                'name': obj['name'],
                'coord': obj['coord'],
                'rise_times': rise_times,
                'set_times': set_times,
                'entry_times': entry_times,
                'exit_times': exit_times,
                'has_dark_event': any(is_dark(Time(time), observer) for times in [rise_times, set_times, entry_times, exit_times] for time in times)
            })

    return results


def get_altaz_at_time(coord, observer, time):
    altaz = coord.transform_to(AltAz(obstime=Time(time), location=observer.location))
    return altaz.alt.deg, altaz.az.deg


def main():
    csv_file = 'astronomical_objects.csv'
    lat, lon = 51.5074, -0.1278  # London
    alt, az, width, height = 60, 250, 160, 60

    timezone = ZoneInfo('Europe/London')
    location = EarthLocation(lat=lat*u.deg, lon=lon*u.deg)
    observer = Observer(location=location, timezone=timezone)

    current_time = datetime.now(timezone)
    start_time = Time(current_time)
    twilight_start, twilight_end = get_twilight_times(observer, start_time, start_time + timedelta(days=1))

    if twilight_end:
        end_time = Time(twilight_end)
    else:
        end_time = start_time + timedelta(days=1)

    objects = read_csv(csv_file)
    visible_objects = calculate_visibility(objects, observer, alt, az, width, height, start_time, end_time, timezone)

    dark_visible_objects = [
        obj for obj in visible_objects
        if any(twilight_start <= time <= end_time for time in obj['entry_times'] + obj['exit_times'])
    ]

    sorted_objects = sorted(dark_visible_objects, key=lambda obj: min(
        min((abs((time - current_time).total_seconds()) for time in obj['entry_times'] if twilight_start <= time <= end_time), default=float('inf')),
        min((abs((time - current_time).total_seconds()) for time in obj['exit_times'] if twilight_start <= time <= end_time), default=float('inf'))
    ))

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
        print(f"{obj['name']}:")
        for times, event_name in [('entry_times', 'Enters'), ('exit_times', 'Exits')]:
            for time in obj[times]:
                if twilight_start <= time <= end_time:
                    event_alt, event_az = get_altaz_at_time(obj['coord'], observer, time)
                    print(f"  {event_name} window: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}*")
                    print(f"    Alt/Az at {event_name.lower()}: {event_alt:.2f}°/{event_az:.2f}°")
        print()


if __name__ == "__main__":
    main()
