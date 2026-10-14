# Sample metadata: GEOME, GBIF, and NCBI

RiboPilot can optionally link each sample to corresponding NCBI, GBIF,
or GEOME records. This information can be shown in the app, flagged
where databases disagree, and used in your FASTA headers for submission
to GenBank.

| Database | What it holds | Identifier you supply |
|----|----|----|
| [GEOME](https://geome-db.org) | Specimen, tissue, and collecting-event records from field projects | **BCID**, such as `ark:/21547/CXu2MBIO1000.1` |
| [GBIF](https://www.gbif.org) | Museum and collection specimen records (“occurrences”) | **gbifID**, such as `6186461308`, or a Smithsonian **EZID** |
| [NCBI BioSample](https://www.ncbi.nlm.nih.gov/biosample) | The biological material behind GenBank and SRA (NCBI’s raw-read archive) data | **BioSample** (`SAMN29555051`) or **SRA** accession (`SRR21844202`) |

RiboPilot never merges values or picks one source over another. Each
source keeps its own fields, and you decide which ones go into your
FASTA headers. Records are not live: if the remote database changes,
fetch again. Fetched data never changes how RiboPilot resolves taxonomy.

**Public records only.** Private GEOME or NCBI projects are not
supported.

Related pages:

- [Linking records across
  databases](https://nmnh-oceandna.github.io/RiboPilot/articles/Metadata-Linking.md):
  fill in the other databases from one ID.
- [NMNH voucher
  submissions](https://nmnh-oceandna.github.io/RiboPilot/articles/NMNH-Vouchers.md):
  NMNH metadata requirements for GenBank records.
- [Metadata field
  reference](https://nmnh-oceandna.github.io/RiboPilot/articles/Metadata-Reference.md):
  every export field, the column names RiboPilot recognizes, and fetch
  messages.

## Quick start

1.  Add a column of IDs to your [mapping
    file](https://nmnh-oceandna.github.io/RiboPilot/articles/Your-Own-Project.html#the-mapping-file).
2.  Name it in
    [`new_project()`](https://nmnh-oceandna.github.io/RiboPilot/reference/new_project.md)
    (for example `mapping_gbif = "Occurrence"`). RiboPilot fetches the
    records during setup.
3.  In the app, check the **Metadata** column: one logo per fetched
    source.
4.  Click the **Metadata** button above a table and tick **Export** on
    the fields you want, such as `gbif_lat_lon` and
    `gbif_collection_date`.
5.  In **Export Data**, click those fields to add them to the FASTA
    header template, then export.

## Common tasks

| Task | In R | In the app |
|----|----|----|
| Link IDs at setup | `new_project(mapping_geome =, mapping_gbif =, mapping_ncbi =)` |  |
| Add or change IDs | [`update_sample_metadata()`](https://nmnh-oceandna.github.io/RiboPilot/reference/update_sample_metadata.md), or `fetch_*(ids =, bcids = / gbifs = / ncbi_ids =)` | Click a sample’s Metadata icon, edit the ID, **Fetch** |
| Refresh records | [`fetch_geome()`](https://nmnh-oceandna.github.io/RiboPilot/reference/fetch_geome.md), [`fetch_gbif()`](https://nmnh-oceandna.github.io/RiboPilot/reference/fetch_gbif.md), [`fetch_ncbi()`](https://nmnh-oceandna.github.io/RiboPilot/reference/fetch_ncbi.md) | **Refresh all** in the viewer |
| Fill in other databases | `link_sources = TRUE` | **Follow links between databases** in the viewer |
| Pick which mapping-file columns are compared | [`set_metadata_columns()`](https://nmnh-oceandna.github.io/RiboPilot/reference/set_metadata_columns.md) | **Mapping file columns…** on the Compare tab |
| Use fields in FASTA headers |  | **Metadata** button, tick **Export** |
| Remove fetched data | [`remove_metadata()`](https://nmnh-oceandna.github.io/RiboPilot/reference/remove_metadata.md) | **Remove fetched data…** in the viewer |

## Finding the identifiers

- **GEOME BCID.** Open a record on [geome-db.org](https://geome-db.org)
  and copy its identifier, which starts with `ark:/`. Use the most
  specific record that matches your sequenced material, usually the
  tissue: RiboPilot also fetches the sample, collecting event,
  expedition, and project above it.
- **gbifID.** Open the occurrence on [gbif.org](https://www.gbif.org)
  and copy the number at the end of its address
  (`https://www.gbif.org/occurrence/6186461308`).
- **Smithsonian EZID.** NMNH specimens can use their EZID, a permanent
  identifier starting with `ark:/65665/`, in place of a gbifID.
  RiboPilot finds the matching GBIF occurrence, so the specimen must be
  on GBIF. EZIDs are more stable than gbifIDs. Links from `n2t.net` or
  `collections.nmnh.si.edu` work, with or without dashes.
- **BioSample or SRA accession.** Copy the BioSample accession (`SAMN`,
  `SAMEA`, or `SAMD`) from [NCBI
  BioSample](https://www.ncbi.nlm.nih.gov/biosample), or use the SRA
  accession of your reads (an `SRR`, `ERR`, or `DRR` run, an `SRX`
  experiment, or an `SRS` sample). For an SRA accession, RiboPilot finds
  its BioSample and also stores the run, experiment, study, platform,
  and library.

Full links work everywhere: RiboPilot trims them to the bare ID.

**gbifIDs can change.** GBIF may give an occurrence a new gbifID when
its dataset is republished. If an ID stops working, RiboPilot keeps the
data it fetched before and marks the fetch as failed. Look the specimen
up again on gbif.org and paste the new ID.

## Adding identifiers

### At project setup

Add a column for any or all databases to your mapping file (see [Your
own
project](https://nmnh-oceandna.github.io/RiboPilot/articles/Your-Own-Project.html#the-mapping-file)).
Columns named `GEOME_BCID`, `GBIF_ID`, and `BioSample` (any upper or
lower case) are found automatically. For other names, pass
`mapping_geome`, `mapping_gbif`, and `mapping_ncbi`:

    ID,Taxon,R1,R2,Tissue_BCID,Occurrence,SRA
    Fish_01,Psenes pellucidus,Fish_01_R1.fastq.gz,Fish_01_R2.fastq.gz,ark:/21547/CXu2MBIO1000.1,,
    Fish_02,Notemigonus crysoleucas,Fish_02_R1.fastq.gz,Fish_02_R2.fastq.gz,,6186461308,
    Fish_03,Conger oceanicus,Fish_03_R1.fastq.gz,Fish_03_R2.fastq.gz,,,SRR21844202

``` r

RiboPilot::new_project(
  path          = "~/my_project",
  mapping_fn    = "~/my_samples.csv",
  data_path     = "~/my_reads",
  mapping_geome = "Tissue_BCID",
  mapping_gbif  = "Occurrence",
  mapping_ncbi  = "SRA"
)
```

- Blank cells are fine: those samples have no link.
- Fetching needs internet access on the machine running R or the app.
  Offline, pass `fetch_geome = FALSE`, `fetch_gbif = FALSE`, and
  `fetch_ncbi = FALSE`, and fetch later. Pipeline jobs never contact
  these databases.

[`new_test_project()`](https://nmnh-oceandna.github.io/RiboPilot/reference/new_test_project.md)
sets up the test project with linking on. Its `fetch_metadata` argument
turns the NCBI fetch off.

### Later

In the app, click the icon in the Metadata column, type or paste the ID,
and click **Fetch**.

In R,
[`add_samples()`](https://nmnh-oceandna.github.io/RiboPilot/reference/add_samples.md)
and
[`update_sample_metadata()`](https://nmnh-oceandna.github.io/RiboPilot/reference/update_sample_metadata.md)
take the same `mapping_*`, `fetch_*`, and `link_sources` arguments. With
[`update_sample_metadata()`](https://nmnh-oceandna.github.io/RiboPilot/reference/update_sample_metadata.md),
samples whose ID changed or was added are fetched again, and a blank
cell removes that sample’s link and data for that source only.

To fetch again without touching the mapping file:

``` r

RiboPilot::fetch_gbif("~/my_project")                   # every gbifID
RiboPilot::fetch_gbif("~/my_project", ids = "Fish_02")  # one sample
RiboPilot::fetch_gbif("~/my_project", ids = c("Fish_01", "Fish_02"),
                      gbifs = c("2336663130", "6186461308"))  # set new IDs, then fetch
RiboPilot::fetch_ncbi("~/my_project", ids = "Fish_03", ncbi_ids = "SAMN29555051")
RiboPilot::fetch_geome("~/my_project", ids = "Fish_01",
                       bcids = "ark:/21547/CXu2MBIO1000.1")
```

A blank ID removes the link. When some fetches fail, the rest still
finish and a warning message lists each failed sample with the reason. A
failed refresh keeps the data fetched before; changing to a different ID
clears the old record right away. See the [fetch
messages](https://nmnh-oceandna.github.io/RiboPilot/articles/Metadata-Reference.html#fetch-messages)
table for additional details.

NCBI allows three requests per second, or ten with a free [NCBI API
key](https://nmnh-oceandna.github.io/RiboPilot/articles/Installation.html#optional-an-ncbi-api-key).
Metadata fetches use the project’s `ncbi_api_key` or the
`NCBI_API_KEY`/`ENTREZ_KEY` environment variables.

## Viewing records in the app

### The Metadata column

The Taxonomy, Assemble, Annotation Review, and Export tables have a
**Metadata** column, with one logo per source the sample has an ID for:

| Icon | Meaning |
|----|----|
| “+” | No IDs. Click it to add one. |
| GEOME “G”, GBIF leaf, NCBI helix | The record was fetched. |
| Faded logo | The ID is set but not fetched yet. |
| Faded logo with a warning triangle | The fetch failed. Hover to see why. |
| Orange flag | The sources disagree on at least one item (see [Comparing sources](#comparing-sources)). |

Hovering over an icon lists each source with its status, then any items
in conflict.

### Showing metadata in the tables

The **Metadata** button above each table opens the **Metadata fields**
window. It lists every field with a value for at least one sample: your
mapping-file columns, then GEOME, GBIF, and NCBI fields, with an example
value and a sample count (for example `12/20`). Each row has two ticks:

- **Show** adds the field as a column in the Taxonomy, Assemble,
  Annotation Review, and Export tables.
- **Export** makes the field usable in FASTA header templates and writes
  it to the sample summary CSV.

Mapping-file columns start with Show ticked and are always available at
export. Click **Save** to keep your choices in the project.

### The sample metadata viewer

Click a sample’s Metadata icon to open the viewer. It has a tab for each
database and a **Compare** tab. Here you can:

- edit an ID and click **Fetch** (clear the box to remove the link),
- click **Refresh all** to fetch every sample again,
- filter fields by name, copy any value, or click **Export CSV** to save
  all of a sample’s metadata.

The **Sample** list on the left marks samples that failed to fetch or
have a conflict. On the **NCBI** tab, records are shown as SRA,
BioSample, and BioProject cards, each linked to its NCBI page.

### Comparing sources

The **Compare** tab puts the values from your mapping file, GEOME, GBIF,
and NCBI side by side for coordinates, collection date, country,
locality, voucher, collector, sex, life stage, and taxon.

| Status | Meaning |
|----|----|
| agree | Every source with a value agrees. Coordinates agree within 0.01 degree. |
| conflict | A real disagreement. |
| not checked | Free text, or a value that could not be read. |
| single | Only one source has a value. |
| \- | No source has a value. |

Blanks never conflict, and placeholders such as `missing` or
`not collected` count as blank. RiboPilot finds your mapping-file
columns by [common
names](https://nmnh-oceandna.github.io/RiboPilot/articles/Metadata-Reference.html#mapping-file-columns-used-by-compare).
If yours differ, pick them under **Mapping file columns…** at the bottom
of the Compare tab, or in R:

``` r

RiboPilot::set_metadata_columns("~/my_project", country = "Country_Name",
                                coordinates = c("Lat_DD", "Long_DD"))
RiboPilot::set_metadata_columns("~/my_project", country = NA)  # back to automatic
RiboPilot::set_metadata_columns("~/my_project", sex = "")      # do not compare sex
```

## Using metadata fields at export

A fetched value reaches your FASTA headers only after you tick
**Export** for it in the **Metadata** button (see [Showing metadata in
the tables](#showing-metadata-in-the-tables)). Each source offers
**GenBank-ready** fields, already formatted as GenBank source modifiers:
the `[name=value]` tags in a FASTA header, such as
`[lat_lon=17.48 S 149.90 W]`. The other fields are the raw values from
each record. See [GenBank
submission](https://nmnh-oceandna.github.io/RiboPilot/articles/GenBank-Submission.md)
for header templates and the modifiers GenBank requires.

In the **Export Data** window, the **Available metadata** panel lists
every field you can use. Click one to insert it into the FASTA header
box. A GenBank-ready field inserts the whole modifier, such as
`[lat_lon={gbif_lat_lon}]`. **Choose fields** next to a source opens the
Metadata fields window so you can enable more fields for export.

A template that takes every specimen value from GBIF:

    {seqid} [organism={organism}] [moltype=genomic DNA] [lat_lon={gbif_lat_lon}] [collection_date={gbif_collection_date}] [geo_loc_name={gbif_geo_loc_name}] [specimen_voucher={gbif_specimen_voucher}] {organism} {rrna_definition}

- **Take a record’s values from one source where you can.** Mixing
  sources can combine coordinates and dates from different collecting
  events. If you cite a BioSample with `[BioSample=]`, GenBank expects
  the header values to agree with it.
- **Empty values are left out.** A modifier with no value for a sample
  is dropped from that sample’s header instead of written empty, which
  GenBank rejects. The check under each header box names the records
  with empty modifiers, and warns when the same modifier appears twice.
  GenBank requires `collection_date` and `geo_loc_name`, so fill them in
  or use a missing-value term such as `missing` (see [Required source
  modifiers](https://nmnh-oceandna.github.io/RiboPilot/articles/GenBank-Submission.html#required-source-modifiers)).
- **Check GBIF values.** GBIF’s `country` is its own English name
  (`United States of America`, where GenBank expects `USA`), and GBIF
  may have adjusted or generalized coordinates. See the [GBIF
  notes](https://nmnh-oceandna.github.io/RiboPilot/articles/Metadata-Reference.html#gbif-notes).
- **Your mapping file wins on a name clash.** If your mapping file has a
  column with the same name as a ticked field (for example
  `geome_lat_lon`), its value is used in headers and the summary.

When you click **Export**, RiboPilot checks the items your templates
use. If any sample in the export group has a **conflict** on one of
them, the **Sample metadata disagrees** warning lists each sample, item,
and value. Click **Export anyway** or **Cancel**.

NMNH specimens need two extra modifiers; see [NMNH voucher
submissions](https://nmnh-oceandna.github.io/RiboPilot/articles/NMNH-Vouchers.md).

## Removing fetched data

**Remove fetched data…** under the sample list in the viewer, or
[`remove_metadata()`](https://nmnh-oceandna.github.io/RiboPilot/reference/remove_metadata.md),
deletes fetched records for one sample or all:

``` r

RiboPilot::remove_metadata("~/my_project")                    # everything fetched
RiboPilot::remove_metadata("~/my_project", sources = "GBIF")  # one source
RiboPilot::remove_metadata("~/my_project", sources = c("GEOME", "GBIF"),
                           ids = "Fish_02")
```

This removes the stored records, their fetch status, and IDs that
linking filled in. Your mapping-file columns, the IDs you supplied, and
your field ticks are kept, so **Refresh all** or `fetch_*()` brings the
data back.
