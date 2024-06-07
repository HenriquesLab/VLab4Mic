<!-- markdownlint-disable -->

<a href="../supra_molecular_simulator/download.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `download`




**Global Variables**
---------------
- **headers**

---

<a href="../supra_molecular_simulator/download.py#L9"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `download_suggested_structures`

```python
download_suggested_structures(data_path: str = 'data') → None
```

Downloads PDB *.cif.gz files for the suggested structures in the given path. 



**Args:**
 
 - <b>`data_path`</b>:  The path to the data directory. Defaults to "data". 



**Returns:**
 None 


---

<a href="../supra_molecular_simulator/download.py#L59"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `download_file`

```python
download_file(url: str, fname: str, chunk_size: int = 1024) → None
```

Downloads a file from the given URL and saves it to the specified file name. 



**Args:**
 
 - <b>`url`</b>:  The URL of the file to download. 
 - <b>`fname`</b>:  The file name to save the downloaded file as. 
 - <b>`chunk_size`</b>:  The size of the chunks to download the file in. Defaults to 1024. 



**Returns:**
 None 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
