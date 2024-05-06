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

        if not os.path.exists(self.outdir):
            os.system(f"mkdir -p {self.outdir}")

    def add_rubbish(self, file):
        self._rubbish.append(os.path.abspath(file))
    
    def add_log_record(self, record):
        self._log_records.append(record)

    def clean(self):
        if len(self._rubbish) == 0:
            return
        else:
            cmds = [f"rm -rf {f}" for f in self._rubbish]
            for cmd in cmds:
                os.system(cmd)
    
    def get_pct(self, val, total):
        return f"{round(100*val/total, 2)} %"

    @staticmethod
    def check_mkdir(dir_name):
        """if dir_name is not exist, make one"""
        if not os.path.exists(dir_name):
            os.system(f"mkdir -p {dir_name}")

    def write_log(self):
        with open(self.log_file, "w") as fd:
            s = "\n".join(self._log_records)
            fd.write(s)
