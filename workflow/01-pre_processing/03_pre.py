### Imports ### 
from pathlib import Path
import pandas as pd

### Constants ### 
base_dir = Path().cwd()
mos_dir = base_dir/"mosdepth"


### Main ###
files = [ i.resolve() for i in mos_dir.glob("*.bed.gz")]
files_done = set()
main_df = pd.DataFrame(columns=['chr', 'start', 'end'])
df = pd.DataFrame(columns=['chr', 'start', 'end'])

for i,j in enumerate(files):
    if main_df.empty: 
        main_df = pd.read_csv(files[i],sep="\t",header=None,names=['chr','start','end',f'depth_{i}'],compression='gzip')
    elif i not in files_done: 
        df = pd.read_csv(files[i],sep="\t",header=None,names=['chr','start','end','mean_depth'],compression='gzip')
    
    if not df.empty:
        df = df.rename(columns={'mean_depth':f'depth_{i}'})
        main_df = pd.merge(main_df,df, on=['chr','start','end'],how = 'left')
        files_done.add(i)

main_df.to_csv(mos_dir/"wide_df.csv",index=False)

