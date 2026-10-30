"""Copy a consistent completed QR block while a bounded diagnostic is running."""
import argparse
import datetime
import hashlib
import io
import json
from pathlib import Path
import numpy as np


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Choose a new checkpoint path')
    metadata=(args.input/'analysis.json').read_bytes()
    history=(args.input/'history.npz').read_bytes()
    record=json.loads(metadata)
    data=np.load(io.BytesIO(history))
    observed=record['observed_tau'];elapsed=record['elapsed_tau']
    d=record['quotient_dimension']
    checks=dict(complete_block=elapsed==200.+observed and observed>0 and observed%500==0,
        same_time=float(data['time'][-1])==elapsed,
        state_endpoint=np.array_equal(data['states'][-1],data['final_state']),
        leading_log=data['log_sums'][0]/observed==record['quotient_leading_ftle'],
        contrast_logs=np.array_equal(data['log_sums'][1:].reshape(-1,3)/observed,np.asarray(record['contrast_ftle'])),
        finite_states=bool(np.isfinite(data['states']).all()),
        tangent_shapes=data['final_quotient_tangent'].shape==(d,) and data['final_contrast_bases'].shape==(len(record['contrast_ftle']),3,3))
    if not all(checks.values()):
        raise ArithmeticError(checks)
    args.output.mkdir(parents=True)
    (args.output/'analysis.json').write_bytes(metadata)
    (args.output/'history.npz').write_bytes(history)
    (args.output/'source.py').write_bytes(Path(__file__).read_bytes())
    result=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),status='passed',checks=checks,
        source_folder=str(args.input),observed_tau=observed,elapsed_tau=elapsed,
        hashes={name:hashlib.sha256((args.output/name).read_bytes()).hexdigest() for name in ['analysis.json','history.npz','source.py']},
        purpose='Freeze an exact complete QR-block boundary. A later wall-budget checkpoint may contain partial-block state without matching accumulated logs and must not be resumed as a complete block.')
    (args.output/'checkpoint.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(result['status'],observed,checks)


if __name__=='__main__':
    main()
