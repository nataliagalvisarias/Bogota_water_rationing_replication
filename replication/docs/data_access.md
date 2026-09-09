# Obtaining the restricted incident data

## What is restricted

Incident-level records from the Bogotá 123 Emergency Line, held by the
**Secretaría Distrital de Seguridad, Convivencia y Justicia** of the Bogotá District
Government. The author obtained these under a data-access request and is not permitted
to redistribute them, including to referees or to a journal data editor.

## How to request an identical extract

Requests are made through the Bogotá District Government's public information channel,
*Bogotá te escucha* — Sistema Distrital para la Gestión de Peticiones Ciudadanas
(https://bogota.gov.co/sdqs/), addressed to the Secretaría Distrital de Seguridad,
Convivencia y Justicia, or by written petition (*derecho de petición*) under Law 1755 of
2015. Colombian law obliges the entity to respond within 15 working days.

### Text of the request that produced the extract used here

> Solicito, en formato de datos abiertos (CSV o similar), el registro de incidentes
> reportados a la Línea de Emergencias 123 de Bogotá, clasificados como delitos contra la
> mujer (feminicidio, delitos sexuales y maltrato a mujer), para el período comprendido
> entre el 1 de enero de 2023 y el 31 de mayo de 2025. Para cada incidente solicito:
> identificador del registro, fecha, hora, tipo de incidente y coordenadas
> georreferenciadas del lugar del hecho. No solicito ningún dato personal identificable
> de las personas involucradas.

### Verifying that a newly obtained extract matches

`docs/checksums.md5` contains checksums for every file shipped in this package. A
replicator who obtains the restricted extract should additionally verify that it matches
the one used here on the following aggregate quantities, all of which are printed by
`code/02_build_panel.R`:

| Quantity | Expected value |
|---|---|
| Domestic-violence reports in the study window, before screens | 138,130 |
| Removed by the adjacency screen | 87,298 (63.2%) |
| Remaining after the adjacency screen | 50,832 |
| Final estimation sample (within 2,000 m) | 20,729 |
| Sexual-violence reports used as the comparison outcome | 2,350 |
| Rationing days in the study period | 325 |

If any of these differs, the extract is not the same vintage and results will not match
exactly. The Secretariat revises historical records, so an extract obtained later may
contain small retrospective corrections.

## Data the Secretariat will not supply

Two inputs required by the extensions in `code/99_extensions_rate_estimand.R` were not
available and were not requested:

1. A **residential population and *estrato* layer** at 50-metre resolution, required for
   the exposure-normalised rate estimand of Section 4.5.2 and for boundary covariate
   balance. Candidate sources are the DANE 2018 census block layer and the Bogotá
   *estrato* cadastral layer; neither is at the resolution the estimand needs without
   areal interpolation, whose own boundary behaviour would need to be defended.
2. **Georeferenced pre-rationing reports for 2023** in the same format, required for the
   placebo-in-time. The 2023 records exist in the extract but with coarser geocoding.
