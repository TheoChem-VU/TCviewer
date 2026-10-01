import json
import os

STYLE_DIRECTORY = os.path.join(os.path.dirname(__file__), 'styles')

class Settings:
	def __init__(self, style_path: str = 'PyOrbb'):
		self.load_style(style_path)
		self.style_path = style_path
		self.is_edited = False

	def load_style(self, style_path: str):
		'''
		Read a style .json file
		'''
		if not os.path.exists(style_path):
			style_path = os.path.join(STYLE_DIRECTORY, style_path)
		if not os.path.exists(style_path):
			style_path = style_path + '.json'
		if not os.path.exists(style_path):
			raise ValueError(f'Could not find path to the style')

		with open(style_path) as sp:
			self.settings = json.loads(sp.read())

	def write_style(self, style_path: str):
		'''
		Write a style .json file
		'''
		style_path_dir = os.path.dirname(style_path)
		if not os.path.exists(style_path_dir):
			style_path = os.path.join(STYLE_DIRECTORY, style_path)

		if style_path.endswith('PyOrbb.json'):
			raise ValueError('Cannot overwrite default PyOrbb.json style')

		if not style_path.endswith('.json'):
			style_path = style_path + '.json'

		with open(style_path, 'w+') as sp:
			sp.write(json.dumps(self.settings))

	def get(self, *keys):
		''' 
		Obtain a value that may be nested in the settings dictionary
		''' 
		d = self.settings
		for key in keys:
			if not isinstance(d, dict):
				raise KeyError(f'Cannot access value at {str(list(keys))}')
			d = d[key]
		return d

	def set(self, *keys):
		''' 
		Set a value that may be nested in the settings dictionary
		'''
		self.is_edited = True
		d = self.settings
		for key in keys[:-2]:
			d = d[key]

		d[keys[-2]] = keys[-1]


s = Settings()
# print(s.settings)
# print(s.get('a', 'b'))
s.set('a', 'heyhey')
print(s.get('a'))

s.write_style('PyOrbb2')

s.load_style('PyOrbb2')
print(s.get('a'))


s.load_style('PyOrbb')
print(s.get('a', 'b'))

