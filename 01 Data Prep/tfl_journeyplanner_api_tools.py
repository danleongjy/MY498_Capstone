import requests
import json
import sqlite3
import traceback
import sys
from math import cos, asin, sqrt
from shapely.geometry import LineString, Point

import tfl_journeyplanner_api_settings as settings

def setup_db(db):
    '''Assumes db is a path to a sqlite3 database.

    Creates db if needed and resets the traveloptions data table.'''
    writedb = sqlite3.connect(db)
    writecursor = writedb.cursor()

    # create ttm table if it doesn't already exist
    writecursor.execute("DROP TABLE IF EXISTS traveloptions")
    writecursor.execute('''CREATE TABLE traveloptions(
                          orig_id TEXT NOT NULL,
                          dest_id TEXT NOT NULL,
                          option_id INT NOT NULL,
                          starttime TEXT NOT NULL,
                          traveltime REAL NOT NULL,
                          fare REAL,
                          transfers INT NOT NULL,
                          legs INT,
                          dist REAL,
                          first_zoned_station TEXT,
                          last_zoned_station TEXT,
                          PRIMARY KEY (orig_id, dest_id, option_id))''')
    writedb.commit()
    writedb.close()

def GenerateODPairs(orig_coords_list, dest_coords_list):
    '''Assumes orig_coords_list and dest_coords_list are lists of tuples (id,x,y) that
    specify the latlong coordinates of each origin and destination respectively.

    Returns a list of all pairs of origin and destination coordinates as long as
    the origin and destination coordinates are different.'''
    odpairs = []
    for orig in orig_coords_list:
        for dest in dest_coords_list:
            if orig == dest:
                continue
            else:
                odpairs.append((orig, dest))
    return odpairs

def call_api(odpair, stoplist, batchsettings):
    '''Assumes odpair is a tuple of two latlong coordinate tuples, stoplist
    is a dictionary with NAPTANs as keys and [x,y] coordinates as values and
    batchsettings is a dictionary of settings for the API call.

    Composes a URI for the TfL Journey Planner API and executes the call.'''
    orig_coord = odpair[0]
    dest_coord = odpair[1]

    # assemble API call
    call = ('https://api.tfl.gov.uk/Journey/JourneyResults/'
            + orig_coord[2] #TfL API expects latlongs, not longlats
            + ','
            + orig_coord[1]
            + '/to/'
            + dest_coord[2]
            + ','
            + dest_coord[1]
            + '?app_id=' + settings.get_app_id()
            + '&app_key=' + settings.get_app_key()
            + batchsettings['date']
            + batchsettings['time']
            + batchsettings['timeIs']
            + batchsettings['journeyPreference']
            + batchsettings['useMultiModalCall']
           )

    # error handling filename
    # errorfilename = 'error-' + str(odpair[0][0]) + ',' + str(odpair[1][0] + '.txt')

    # execute call
    try:
        result = requests.get(call, timeout = (10, 600))
    # if bad response from requests module
    except Exception as e:
        # with open(errorfilename, 'w') as errorfile:
        #     traceback.print_exc(file = errorfile)
        return [[odpair, None, None, None, False]]

    # if successful
    if result.status_code == 200:
        # load JSON data
        dboutputs = []

        data = json.loads(result.content)
        # for each travel option in the results
        for option in range(len(data['journeys'])):
            # parse the JSON data and add it to outputs
            option_output = parse_json(data['journeys'][option])
            network_output = parse_network(data['journeys'][option], stoplist)
            dboutputs.append([odpair, option, option_output, network_output, True])
        return dboutputs
    # if server-side error
    else:
        # with open(errorfilename, 'w') as errorfile:
        #     errorfile.write(result.text)
        return [[odpair, None, None, None, False]]

def calculate_distance(coord_list):
    '''Assumes coord_list is the output of parse_linestring.

    Haversine function implementation from
    https://stackoverflow.com/questions/27928/calculate-distance-between-two-latitude-longitude-points-haversine-formula.
    Modified to return metres.'''

    distance = 0
    p = 0.017453292519943295 # pi / 180, for conversion from deg to rad
    for i in range(len(coord_list) - 1):
        a = 0.5 - cos((float(coord_list[i + 1][1]) - float(coord_list[i][1])) * p)/2 + cos(float(coord_list[i][1]) * p) * cos(float(coord_list[i + 1][1]) * p) * (1 - cos((float(coord_list[i + 1][0]) - float(coord_list[i][0])) * p)) / 2
        distance += 12742 * asin(sqrt(a))

    return distance * 1000

def parse_linestring(linestring):
    '''Assumes linestring is a lineString from TfL JSON outputs.

    Parses the linestring into a list of (x,y) coordinates.'''
    coord_list = []
    for coord in linestring[2:-2].split('],['):
        coordlat = float(coord.split(', ')[0]) # TfL JSON outputs latlongs
        coordlong = float(coord.split(', ')[1])
        if [coordlong, coordlat] not in coord_list:
            coord_list.append([coordlong, coordlat]) # standardise to longlats
    return coord_list

def parse_json(data):
    '''Assumes data is a TfL Journey Planner API journey option result in JSON format.

    Parses the JSON outputs to produce the data for writing to the database.'''
    outputs = {'starttime': data['startDateTime'],
               'traveltime': data['duration'],
               'fare': None,
               'transfers': 0,
               'legs': 0,
               'dist': 0,
               'first_zoned_station': None,
               'last_zoned_station': None}

    if 'fare' in data.keys():
        outputs['fare'] = data['fare']['totalCost']

    for leg in data['legs']:
        # try to calculate the leg distance from the path linestring
        try:
            leg_dist = calculate_distance(parse_linestring(leg['path']['lineString']))
        # if there is no linestring, calculate the leg distance from the departure and arrival points
        except KeyError:
            leg_dist = calculate_distance(parse_linestring('[' +
                                                          str([leg['departurePoint']['lon'],leg['departurePoint']['lat']]) +
                                                          ',' +
                                                          str([leg['arrivalPoint']['lon'],leg['arrivalPoint']['lat']]) +
                                                          ']'))
        leg_mode = leg['mode']['id']
        leg_time = leg['duration']

        outputs['legs'] +=1
        outputs['dist'] += leg_dist
        # add data on each leg by mode to the outputs dict
        if 'dist_' + leg_mode in outputs.keys():
            outputs['dist_' + leg_mode] += leg_dist
        else:
            outputs['dist_' + leg_mode] = leg_dist

        if 'legs_' + leg_mode in outputs.keys():
            outputs['legs_' + leg_mode] += 1
        else:
            outputs['legs_' + leg_mode] = 1

        if 'time_' + leg_mode in outputs.keys():
            outputs['time_' + leg_mode] += leg_time
        else:
            outputs['time_' + leg_mode] = leg_time

        # if the mode is one of the motorised ones, increment transfers
        if leg_mode in {'tube', 'overground', 'tflrail', 'national-rail', 'dlr', 'tram', 'bus', 'coach', 'river-bus'}:
            outputs['transfers'] += 1

        # if the mode is one of those that uses farezones, try to identify the first and last zoned station naptans
        if leg_mode in {'tube', 'overground', 'tflrail', 'national-rail', 'dlr', 'tram'}:
            if outputs['first_zoned_station'] is None:
                try:
                    outputs['first_zoned_station'] = leg['departurePoint']['naptanId']
                except KeyError:
                    pass
            try:
                outputs['last_zoned_station'] = leg['arrivalPoint']['naptanId']
            except KeyError:
                pass
    # transfers = motorised legs - 1
    outputs['transfers'] -= 1
    if outputs['transfers'] < 0:
        outputs['transfers'] = 0
    return outputs

def parse_network(data, stoplist):
    '''Assumes data is a TfL Journey Planner API journey option result in JSON format, stoplist
    is a dictionary with NAPTANs as keys and [x,y] coordinates as values, and edges and nodes are
    the output of parse_network.

    Parses the JSON data to obtain data on network edges and appends it to edges and nodes.'''

    edges = {}
    nodes = {}

    for leg in data['legs']:
        mode = leg['mode']['id']
        leg_lines = {leg['routeOptions'][i]['name']: leg['routeOptions'][i]['directions'] for i in range(len(leg['routeOptions']))}

        # compile a list of stoppoint coordinates
        stoppoints = []
        # first stop is the departure point. Use DfT stop coords by default (using departure point NAPTAN if it exists)
        # otherwise use TfL coords.
        if ('naptanId' in leg['departurePoint'].keys() and leg['departurePoint']['naptanId'] in stoplist.keys()):
            stoppoints.append(stoplist[leg['departurePoint']['naptanId']])
        else:
            stoppoints.append([leg['departurePoint']['lon'], leg['departurePoint']['lat']])
        # TfL JSON does not include coords for stopPoints, so we obtain the coords from DfT NAPTAN XML.
        # if the stopPoints list has a length more than 0 and the NAPTANs match the set from DfT,
        # append stopPoint coordinates from the DfT XML. Otherwise omit the stopPoint.
        if len(leg['path']['stopPoints']) > 0:
            for stop in leg['path']['stopPoints'][:-1]: # omit the last stopPoint as we will use the coords from arrivalPoint
                if 'id' in stop.keys():
                    if stop['id'] in stoplist.keys():
                        if stoplist[stop['id']] not in stoppoints:
                            stoppoints.append(stoplist[stop['id']])
        # last stop is the arrival point. Use DfT stop coords by default (using departure point NAPTAN if it exists)
        # otherwise use TfL coords.
        if ('naptanId' in leg['arrivalPoint'].keys() and leg['arrivalPoint']['naptanId'] in stoplist.keys()):
            stoppoints.append(stoplist[leg['arrivalPoint']['naptanId']])
        else:
            stoppoints.append([leg['arrivalPoint']['lon'], leg['arrivalPoint']['lat']])

        # parse the path linestring into a list of coordinates
        if 'lineString' in leg['path'].keys():
            path = parse_linestring(leg['path']['lineString'])
        # if there is no linestring, the path is just the list of stoppoints.
        else:
            path = stoppoints

        # compile a list of subpaths
        # a subpath is the part of a path that is between 2 consecutive stops
        # assume all parts of the path must belong to one subpath
        # this procedure is robust to cases where there is more than one stop between
        # a pair of path points
        subpaths = []
        # create a copy of path. in processing, we will whittle down the path segments
        # as we assign them to subpaths.
        remaining_path = [pt for pt in path]
        # we know that the first stop must be before the first path point
        remaining_path.insert(0, stoppoints[0])
        last_stop_idx = 0

        # consider all subsequent stops in order
        for stop_idx in range(1,len(stoppoints) - 1):
            # start constructing the subpath
            subpath = []

            # consider the current stop. calculate the perpendicular dist to all lines connecting two consecutive points
            # in the remaining path.
            stop_prox_to_segment = []
            for point_idx in range(len(remaining_path) - 1):
                stop_prox_to_segment.append(Point(stoppoints[stop_idx][0],stoppoints[stop_idx][1]).distance(LineString([remaining_path[point_idx], remaining_path[point_idx + 1]])))

            # then linear search to find which segment of the remaining path it is closest to
            current_closest_segment = len(stop_prox_to_segment)
            # check if stop_prox_to_segment is empty. this could happen if there are only 2 path points.
            if len(stop_prox_to_segment) > 0:
                current_mindist_to_segment = max(stop_prox_to_segment)
                for segment_idx in range(len(stop_prox_to_segment)):
                    if stop_prox_to_segment[segment_idx] < current_mindist_to_segment:
                        current_mindist_to_segment = stop_prox_to_segment[segment_idx]
                        current_closest_segment = segment_idx

            # then extend the subpath by the parts of the path that are between the previous path point and the path point of
            # the closest segment, and add the current stop. if current_closest_segment = 0, this will just add the previous
            # stop and the current stop to the subpath.
            subpath.extend(remaining_path[:current_closest_segment + 1])
            subpath.append(stoppoints[stop_idx])

            # append the completed subpath to subpaths, move to the next stop and omit the part of the path that is already in
            # a subpath from the next search. add the current stop to the start of the remaining path so it will be included
            # in the next subpath.
            subpaths.append(subpath)
            remaining_path = remaining_path[current_closest_segment + 1:]
            remaining_path.insert(0, stoppoints[stop_idx])
            last_stop_idx = stop_idx

        # for the last subpath, just add the remaining path points and the last stop
        subpath = [pt for pt in remaining_path]
        subpath.append(stoppoints[-1])
        subpaths.append(subpath)

        # for each subpath, calculate the distance covered and write data to edges and nodes
        for subpath in subpaths:
            distance = calculate_distance(subpath)
            o = str(subpath[0][0]) + ',' + str(subpath[0][1])
            d = str(subpath[-1][0]) + ',' + str(subpath[-1][1])
            od = o + '_' + d + '_' + mode
            if distance > 0:
                # add node data
                if o not in nodes.keys():
                    nodes[o] = {'location': subpath[0], 'routes_via_me': 1}
                else:
                    nodes[o]['routes_via_me'] += 1

                if d not in nodes.keys():
                    nodes[d] = {'location': subpath[-1], 'routes_via_me': 1}
                else:
                    nodes[d]['routes_via_me'] += 1

                # add data to edges if it doesn't already exist
                # edges are unique by od and mode - thus there can be multiple edges with the same od but different mode
                if od not in edges.keys():
                    edges[od] = {'origin': subpath[0], 'destination': subpath[-1], 'mode': mode, 'distance': distance, 'path': subpath, 'routes': {}}

                # add new routing data to edge if it is available
                for route in leg_lines.keys():
                    if route in edges[od]['routes'].keys():
                        edges[od]['routes'][route].union(set(leg_lines[route]))
                    else:
                        edges[od]['routes'][route] = set(leg_lines[route])

        # diagnostics if needed
        ''''
        import matplotlib.pyplot as plt

        print('Mode:', mode)
        print('Path:', len(path))
        print('Stops:', len(stoppoints))
        print('Subpaths:')
        subpathsmap = []
        for subpath in subpaths:
            subpathmap = []
            for coord in subpath:
                if coord in stoppoints:
                    subpathmap.append('S' + str(stoppoints.index(coord)))
                elif coord in path:
                    subpathmap.append('P' + str(path.index(coord)))
                else:
                    subpathmap.append('Err')
            subpathsmap.append(subpathmap)
        print(subpathsmap)

        fig1 = plt.figure(figsize = (10,10))
        ax1 = fig1.add_subplot(111)
        plt.plot([pt[0] for pt in path], [pt[1] for pt in path], label = 'Path')
        for pt_idx in range(len(path)):
            plt.text(path[pt_idx][0], path[pt_idx][1], 'P' + str(pt_idx))
        plt.scatter([pt[0] for pt in path], [pt[1] for pt in path], label = 'PathVertex')
        plt.scatter([pt[0] for pt in stoppoints], [pt[1] for pt in stoppoints], marker = 'X', label = 'Stops')
        for subpath in subpaths:
            plt.plot([pt[0] for pt in subpath], [pt[1] for pt in subpath], '--', label = 'Subpaths')
        plt.legend()
        plt.show()
        '''
    return edges, nodes

def write_to_db(db, orig_id, dest_id, option, json_output_dict):
    '''Assumes db is a sqlite3 database, orig_id and dest_id are unique ids
    for the origin and destination coordinates, option is a unique id for
    a travel option between orig_id and dest_id, and json_output_dict is the
    output of the parse_json function.

    Writes data to the travel time options table.'''
    writedb = sqlite3.connect(db)
    writecursor = writedb.cursor()

    # check which columns are currently in the table
    schema = writecursor.execute('PRAGMA table_info(traveloptions)').fetchall()
    cols = {col[1] for col in schema}

    params = (orig_id, dest_id, option, json_output_dict['starttime'], json_output_dict['traveltime'],
              json_output_dict['fare'], json_output_dict['transfers'], json_output_dict['legs'],
              json_output_dict['dist'], json_output_dict['first_zoned_station'],
              json_output_dict['last_zoned_station'])

    # add in a new record with the default fields
    writecursor.execute('''INSERT INTO traveloptions(orig_id, dest_id, option_id, starttime, traveltime,
    fare, transfers, legs, dist, first_zoned_station, last_zoned_station)
    VALUES(?,?,?,?,?,?,?,?,?,?,?)''',params)

    # identify modes in current legs
    modes = []
    for key in json_output_dict.keys():
        if (len(key) > 5) and ((key[:5] in {'legs_', 'dist_', 'time_'})):
            modes.append(key)
    for mode in modes:
        # if the current mode is not in the table columns, add it as a new column
        if mode not in cols:
            writecursor.execute('ALTER TABLE traveloptions ADD "' + mode + '" REAL DEFAULT 0')
            writedb.commit()
            # update the set of table columns
            cols.add(mode)

        # update the current record with data on the current mode
        writecursor.execute('UPDATE traveloptions SET "' + mode + '" = ' + str(json_output_dict[mode]) +
                            " WHERE orig_id = '" + orig_id + "' AND dest_id = '" + dest_id +
                            "' AND option_id = " + str(option))
        writedb.commit()
    writedb.close()

def write_to_network(edges, nodes, edge_node_tuple):
    '''
    Assumes edges and nodes are dictionaries in the same format as the output of parse_network,
    and edge_node_tuple is the output of parse_network.

    Merges the data from edge_node_tuple into edges and nodes.
    '''

    for edge in edge_node_tuple[0].keys():
        if edge in edges.keys():
            for route in edge_node_tuple[0][edge]['routes'].keys():
                if route in edges[edge]['routes'].keys():
                    edges[edge]['routes'][route] = edges[edge]['routes'][route].union(edge_node_tuple[0][edge]['routes'][route])
                else:
                    edges[edge]['routes'][route] = edge_node_tuple[0][edge]['routes'][route]
        else:
            edges[edge] = edge_node_tuple[0][edge]

    for node in edge_node_tuple[1].keys():
        if node in nodes.keys():
            nodes[node]['routes_via_me'] += edge_node_tuple[1][node]['routes_via_me']
        else:
            nodes[node] = edge_node_tuple[1][node]
