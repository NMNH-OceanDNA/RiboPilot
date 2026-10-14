# Metadata field reference

Look-up tables for [Sample
metadata](https://nmnh-oceandna.github.io/RiboPilot/articles/Sample-Metadata.md).

## GenBank-ready fields

Each field is named `<source>_<name>`, for example `gbif_lat_lon`. A
field is empty for a sample unless its required parts exist. For GEOME,
each value comes from the nearest record that has it (tissue, then
sample, then event). GBIF fields come from the occurrence only, never
from its dataset or publisher. NCBI fields pass the BioSample value
through, since BioSample attributes already use GenBank’s modifier
names.

| Name | GEOME (`geome_`) | GBIF (`gbif_`) | NCBI (`ncbi_`) |
|----|----|----|----|
| `lat_lon` | `decimalLatitude`, `decimalLongitude` | `decimalLatitude`, `decimalLongitude` | `lat_lon` (a decimal pair is converted) |
| `collection_date` | `yearCollected`, `monthCollected`, `dayCollected` | `year`, `month`, `day`; else a single-date `eventDate` | `collection_date` |
| `geo_loc_name` | `country`, then `stateProvince` and `locality` | `country`, then `stateProvince` and `locality` | `geo_loc_name` |
| `specimen_voucher` | `genbankSpecimenVoucher`; else `catalogNumber`, with `institutionCode` as a prefix | `institutionCode`, `collectionCode`, `catalogNumber` | `specimen_voucher`; else `genbankSpecimenVoucher` |
| `collected_by` | `collectorList` | `recordedBy` | `collected_by` |
| `identified_by` |  | `identifiedBy` | `identified_by` |
| `tissue_type` | `tissueType` |  |  |
| `sex` | `sex` | `sex` (lowercased) | `sex` (lowercased) |
| `dev_stage` | `lifeStage` | `lifeStage` (lowercased) | `dev_stage` (lowercased) |

Example values: `lat_lon` `17.48260 S 149.89990 W`; `collection_date`
`2006-03-18`, `2009-11`, or `2009`; `geo_loc_name`
`French Polynesia: Moorea`; `specimen_voucher` `USNM:FISH:407263`.

- **Vouchers.** A voucher or catalog number that is a web address or an
  ARK gives an empty voucher. GBIF vouchers drop a leading repeat of the
  institution code (`USNM 419933` becomes `USNM:FISH:419933`). GenBank
  checks the institution and collection codes against [NCBI
  BioCollections](https://www.ncbi.nlm.nih.gov/biocollections); GBIF’s
  collection codes (such as `Fishes`) often differ from NCBI’s.
- **NCBI placeholders** give an empty value: `missing`, `not collected`,
  `not applicable`, `not provided`, `restricted access`, `not recorded`,
  `not available`, `unspecified`, `unknown`, `none`, `na`, `n/a`, and
  `-`, including forms like `missing: control sample`. A sex of
  `unknown` is therefore empty.

NCBI also offers fields that insert the `[BioSample=]` and
`[BioProject=]` header tags, which GenBank uses to link the record to
the BioSample and BioProject:

| Name                  | Built from                        | Example        |
|-----------------------|-----------------------------------|----------------|
| `ncbi_biosample`      | BioSample accession               | `SAMN29555051` |
| `ncbi_bioproject`     | first BioProject of the BioSample | `PRJNA720393`  |
| `ncbi_biosample_sra`  | the SRA record’s `biosample`      | `SAMN29555051` |
| `ncbi_bioproject_sra` | the SRA record’s `bioproject`     | `PRJNA720393`  |

A BioSample can belong to several BioProjects; check that
`ncbi_bioproject` is the one for your submission. The `_sra` fields are
filled only for samples fetched through an SRA accession.

## Raw fields

Every other field in a fetched record is available as
`<source>_<Level>_<field>`, for example `geome_Tissue_tissueType`,
`gbif_Occurrence_occurrenceID`, or `ncbi_BioSample_isolate`.

## GBIF notes

GBIF interprets the records it publishes, so check its values before
submitting:

- `gbif_lat_lon` is empty when GBIF flags the coordinates
  `ZERO_COORDINATE`, `COORDINATE_INVALID`, or `COORDINATE_OUT_OF_RANGE`.
- Other flags keep the value but are shown as badges in the viewer.
  Check values flagged `PRESUMED_SWAPPED_COORDINATE`,
  `PRESUMED_NEGATED_LATITUDE` or `_LONGITUDE`, `COORDINATE_REPROJECTED`,
  `COORDINATE_ROUNDED`, or `COUNTRY_COORDINATE_MISMATCH`.
- Coordinates of sensitive species may be deliberately generalized.
  Check the `dataGeneralizations` and `coordinateUncertaintyInMeters`
  raw fields.
- GBIF’s `country` is its own English name (`United States of America`),
  not always the name GenBank expects (`USA`). Check every source’s
  `geo_loc_name` against GenBank’s [country
  list](https://www.ncbi.nlm.nih.gov/genbank/collab/country/); marine
  samples use ocean or sea names from the same list.
- `identifiedBy` and the GBIF taxon describe the museum specimen, which
  may differ from your identification of the sequenced tissue.

## Mapping-file columns used by Compare

The Compare tab finds your mapping-file columns by name, ignoring case.
The first match wins. Use
[`set_metadata_columns()`](https://nmnh-oceandna.github.io/RiboPilot/reference/set_metadata_columns.md)
or **Mapping file columns…** when yours differ.

| Item | Columns tried |
|----|----|
| coordinates | `lat_lon`; or a latitude column (`lat`, `latitude`, `decimalLatitude`) plus a longitude column (`lon`, `long`, `longitude`, `decimalLongitude`) |
| collection_date | `collection_date`, `date`, `eventDate` |
| country | `country`, `geo_loc_name` (the part before `:`) |
| locality | `locality` |
| voucher | `specimen_voucher`, `voucher`, `catalogNumber` |
| collector | `collected_by`, `collector`, `recordedBy` |
| sex | `sex` |
| dev_stage | `dev_stage`, `life_stage`, `lifeStage` |
| taxon | `Taxon` (always) |

NCBI’s country and locality come from the two parts of `geo_loc_name`,
and its taxon from the BioSample organism.

## Voucher search fields

When
[linking](https://nmnh-oceandna.github.io/RiboPilot/articles/Metadata-Linking.html#matching-a-voucher-to-a-gbif-record)
searches GBIF for a voucher, it reads:

| Record | Fields |
|----|----|
| NCBI BioSample | `specimen_voucher`, `genbankSpecimenVoucher`, `materialSampleID`, `voucherCatalogNumber`, `bio_material` |
| GEOME | `genbankSpecimenVoucher`, `materialSampleID`, `voucherCatalogNumber`, `otherCatalogNumbers`, `catalogNumber`, or a bare `voucherCatalogNumber` with the `institutionID` |

## Fetch messages

Shown in the viewer, when hovering a faded logo, or in the warning after
a fetch.

| Message | Meaning |
|----|----|
| `BCID not found in GEOME` | No GEOME record has that ID. A typo in a well-formed BCID shows up this way. |
| `record is private or needs a GEOME login` | The record belongs to a private GEOME project, which RiboPilot cannot read. |
| `BCID not recognized by GEOME` | GEOME rejected the ID as malformed. |
| `'x' is not a GEOME BCID (...)` | The value does not look like a BCID. |
| `GBIF occurrence not found (IDs can change when a dataset is republished)` | No GBIF occurrence has that gbifID. Look the specimen up again on gbif.org. |
| `no GBIF occurrence found for NMNH EZID ...` | GBIF has no specimen with that EZID, or has not indexed it yet. Try again later, or use the gbifID. |
| `'x' is not a GBIF occurrence ID (...)` | The value is not a gbifID, gbif.org link, or NMNH EZID. |
| `BioSample ... not found`, `SRA accession ... not found` | NCBI has no record with that accession. |
| `SRA accession ... has no linked BioSample` | Supply the BioSample accession instead. |
| `'x' is not a BioSample or SRA accession (...)` | The value is not a recognized NCBI accession. |
| `could not reach GEOME / GBIF / NCBI (...)` | No connection. Try again later with `fetch_*()` or **Refresh all**. |
| `GEOME / GBIF / NCBI returned HTTP ...` | Any other server error, often a busy server. Try again later; for NCBI, set an [API key](https://nmnh-oceandna.github.io/RiboPilot/articles/Installation.html#optional-an-ncbi-api-key). |
| `N possible GBIF matches for ...; not linked` | Linking found several GBIF specimens for a voucher. Paste the right gbifID. |
| `too many GBIF matches for ...; not linked` | The voucher matched more than 50 GBIF records. Paste the right gbifID. |
| `... record links to ...; kept your ID ...` | A fetched record names a different ID than yours. Check which is right. |
