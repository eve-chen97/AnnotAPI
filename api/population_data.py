#  To accept a list of dictionaries
# Currently not used

import csv
import json

class PopulationData:
    def __init__(self, data=None):
        """
        Initializes a new instance of the PopulationData class.
        
        :param data: Data containing population data which could be a list of dictionaries.
        """
        self.data = data or []
        self.marker_names = []
        if self.data:
            self._load_data()

    def _load_data(self):
        """
        Private method to load data from the given list of dictionaries.
        """
        # Assuming all dictionaries have the same keys
        self.marker_names = list(self.data[0].keys())

    def from_csv(self, file_path):
        """
        Method to load data from a CSV file.
        """
        with open(file_path, mode='r', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            self.marker_names = reader.fieldnames
            self.data = [row for row in reader]

    def from_json(self, file_path):
        """
        Method to load data from a JSON file.
        """
        with open(file_path, 'r') as jsonfile:
            self.data = json.load(jsonfile)
            self._load_data()

    def get_marker_data(self, marker_name):
        """
        Gets the data for a specific marker.
        
        :param marker_name: Name of the marker to get data for.
        :return: List of values for the specified marker.
        """
        return [float(item[marker_name]) for item in self.data if marker_name in item]
