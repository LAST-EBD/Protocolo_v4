import os
import arrow
import landsatxplore.api
from landsatxplore.earthexplorer import EarthExplorer

api = landsatxplore.api.API(os.environ.get('usgs_user'), os.environ.get('usgs_password'))
ee = EarthExplorer(os.environ.get('usgs_user'), os.environ.get('usgs_password'))

today = arrow.now().format('YYYY-MM-DD')
month = arrow.now().shift(days=-45).format('YYYY-MM-DD')

scenes_oli = api.search(
        dataset='LANDSAT_8_C1',
        latitude=37.47944444,
        longitude=-6.25861111,
        start_date=month,
        end_date=today,
        max_cloud_cover=100,
        max_results=1000)

scenes_etm = api.search(
        dataset='LANDSAT_ETM_C1',
        latitude=37.47944444,
        longitude=-6.25861111,
        start_date=month,
        end_date=today,
        max_cloud_cover=100,
        max_results=1000)

scenes_tm = api.search(
         dataset='LANDSAT_TM_C1',
         latitude=37.47944444,
         longitude=-6.25861111,
         start_date=month,
         end_date=today,
         max_cloud_cover=100,
         max_results=1000)


for scene in  scenes_etm:
        sc = scene['displayId']
        if sc.split('_')[2] == '202034' and sc.split('_')[-1] == 'T1':
                #print('Downloading', sc)
                ee.download(scene_id=sc, output_dir='/root/protocolo/ori')


ee.logout()
api.logout()

