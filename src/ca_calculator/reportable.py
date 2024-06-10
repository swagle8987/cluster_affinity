import pandas as pd


class Reportable:
    """
    Records a report if logging value is set to 1
    """

    def __init__(self, logging_level):
        self.logging_level = logging_level
        if self.logging_level == 1:
            self.log = []
        else:
            pass

    def log_event(self, eventID, data):
        if self.logging_level:
            self.log.append({"eventID": eventID, "data": data})
        else:
            pass

    def get_event_dataframe(self, eventID):
        if self.logging_level == 1:
            dl = []
            for i in self.log:
                if i["eventID"] == eventID:
                    dl.append(i["data"])
            df = pd.DataFrame(dl)
            return df
        else:
            pass
