import os

class MinAna:
    def __init__(
            self,
            sample,
            outdir) -> None:
        self.sample = sample
        self.outdir = outdir
        
        self.log_file = f"{self.outdir}/{self.sample}.log"
        self._rubbish = []
        self._log_records = []

    def add_rubbish(self, file):
        self._rubbish.append(file)
    
    def add_log_record(self, record):
        self._log_records.append(record)

    def clean(self):
        if len(self._rubbish) == 0:
            return
        else:
            cmds = [f"rm -rf {f}" for f in self._rubbish]
            for cmd in cmds:
                os.system(cmd)
    
    def write_log(self):
        with open(self.log_file, "w") as fd:
            s = "\n".join(self._log_records)
            fd.write(s)
