# Convert Python scripts
find . -iname "*.py" -type f -exec sh -c 'echo "Converting {}" && mkdir -p ../notebooks/$(dirname {}) && python jupyter_converter.py {} ../notebooks/$(dirname {})/$(basename {} ".py").ipynb' \;

# Copy mesh files
find . -iname "*.xml" -type f -exec sh -c 'echo "Copying {}" && cp {} ../notebooks/$(dirname {})' \;
find . -iname "*.xml.gz" -type f -exec sh -c 'echo "Copying {}" && cp {} ../notebooks/$(dirname {})' \;

# Copy over required Python modules manually
cp 02_static_linear_pdes/kul/plotslopes.py ../notebooks/02_static_linear_pdes/kul
