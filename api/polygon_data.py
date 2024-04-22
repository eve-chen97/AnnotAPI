# To create a JSON payload from an API request.

import json

class PolygonData:
    def __init__(self, vertices=None):
        """
        Initializes a new instance of the PolygonData class.
        :param vertices: A list of vertex data to initialize the object with. Each vertex is expected to be a point that contributes to forming a closed polygon.
        """
        if vertices is None:
            vertices = []
        self.vertices = vertices

    @classmethod
    def from_json(cls, json_data):
        """
        Creates a new instance of PolygonData from JSON data.
        
        :param json_data: A JSON object containing vertices data for the polygon.
        :return: A new instance of PolygonData initialized with data from the JSON object.
        """
        return cls(json_data)

    @classmethod
    def from_json_file(cls, file_path):
        """
        Creates a new instance of PolygonData from a .json file.
        
        :param file_path: The path to the .json file containing vertices data for the polygon.
        :return: A new instance of PolygonData initialized with data from the .json file.
        """
        with open(file_path, 'r') as file:
            vertices = json.load(file)
        return cls(vertices)

    def get_vertex(self, index):
        """
        Retrieves the vertex data at the specified index from the polygon vertices list.
        
        :param index: The index of the vertex to retrieve.
        :return: The data of the vertex at the specified index.
        """
        if 0 <= index < len(self.vertices):
            return self.vertices[index]
        else:
            print("Index out of range")
            return None

    def get_all_vertices(self):
        """
        Retrieves all vertices data, which can be used to form the polygon.
        
        :return: A list of all vertices.
        """
        return self.vertices

    
# Usage within Flask app, assuming JSON data is provided in the request
# json_data would be the parsed JSON data from a client's request
# polygon_data_instance = PolygonData.from_json(json_data)