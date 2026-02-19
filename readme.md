### Disclaimer

```
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
DEALINGS IN THE SOFTWARE.
```
# hgene

**h**erpesvirus **gene** variant analysis from Nanopore sequencing data.

### Usage

`hgene -v <virus> [-c <cpu>] <prefix|prefix.fastq>`

Arguments:

- `v- <virus>` — Virus reference key (e.g. HHV1, HHV2)
- `-c <cpu>` — Number of threads (default: nproc)
- `<input>` — Prefix or an uncompressed .fastq file

Notes:

- .fq files are not supported
- Compressed files (.gz) are not supported

### Acknowledgments

Julien Prados (@pradosj, Bioinformatics Support Platform) — for his contribution.
