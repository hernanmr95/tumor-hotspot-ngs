# Alcance y límites de este proyecto

Este documento existe como artefacto propio (no como una sección más del
README) porque la distinción que describe es la parte más importante del
repositorio: qué es esta herramienta y, sobre todo, qué **no** es.

## Qué es

Una demostración técnica de un componente muy concreto de un pipeline de
tumor profiling: el cruce de variantes ya anotadas contra un panel de
hotspots conocidos, con generación de un informe trazable. Los hotspots
incluidos (BRAF p.Val600Glu, EGFR p.Leu858Arg, KRAS p.Gly12Asp, etc.) son
marcadores clínicamente establecidos y ampliamente documentados en la
literatura de oncología molecular — no son inventados.

## Qué NO es

**No es un dispositivo de diagnóstico clínico**, ni pretende serlo. Para
que una herramienta como esta pudiera usarse en la práctica clínica real
necesitaría, como mínimo:

- **Validación analítica formal**: sensibilidad, especificidad, límite de
  detección y reproducibilidad inter-ensayo, documentadas siguiendo guías
  como las del College of American Pathologists (CAP) o equivalentes.
- **Certificación regulatoria**: CE-IVD en la Unión Europea, aprobación
  FDA en EE.UU., o funcionar como *laboratory-developed test* bajo un
  marco de calidad como CLIA / ISO 15189.
- **Panel de referencia con cobertura clínica real**: cientos o miles de
  variantes (no 11), incluyendo indels, fusiones génicas, variantes de
  número de copia y splicing — y mecanismos para distinguir variantes
  somáticas de germinales.
- **Curación experta continua**: siguiendo marcos de tiering de evidencia
  como las guías conjuntas AMP/ASCO/CAP, actualizados por comités, no por
  un CSV estático.
- **Trazabilidad y control de versiones del panel** con firma de un
  responsable clínico/científico ante cada actualización.

## Decisiones de diseño que reflejan estos límites

- El matching es por gen + cambio de proteína, explícitamente documentado
  como una simplificación (ver docstring de `hotspot_matcher.py`).
- Las variantes con VAF baja no se descartan automáticamente: se reportan
  aparte para revisión humana, porque esa es una decisión que no debería
  tomar el software solo.
- Todas las variantes fuera de panel se listan igualmente en el informe
  (trazabilidad completa), en vez de solo mostrar los "positivos".
- El informe generado incluye un aviso explícito en la cabecera.

## Uso previsto

Portfolio técnico, aprendizaje propio, y como punto de partida conceptual
para conversaciones sobre cómo se diseñaría el componente equivalente en
un entorno con las licencias y la infraestructura de calidad adecuadas.
