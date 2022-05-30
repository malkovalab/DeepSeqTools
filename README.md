# DeepSeqTools
Code used for the processing of the deep sequencing data of lys-2-InsH region. Made to work on UIOWA Argon cluster. The flowchart explaining the processing pipeline and environment is shown below:
![](https://github.com/malkovalab/DeepSeqTools/blob/main/deepSeqPipeline.png)

## Dependencies

This `python3` project relies on the following list of packages and dependencies:
- `matplotlib`
- `pysam`
- `pandas`
- `seaborn`
- `regex`
- `numpy`
- `joblib`

To install the required packages run this command:

```bash
pip install -r requirements.txt
```
The `python` code used for the reads processing as well as scripts used for reads pre-processing (made for uiowa ARGON SGE compute cluster) is placed inside of the code directory.
