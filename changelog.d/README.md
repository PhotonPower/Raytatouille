# Changelog-Fragmente

Jeder PR legt hier eine eigene Datei an, statt `CHANGELOG.md` zu ändern. So stoßen parallele PRs
nicht mehr im CHANGELOG zusammen.

- Name: `<issue>.<art>.md`, bei mehreren Einträgen derselben Art `<issue>.<art>.<n>.md`.
- `<art>`: `added` (Hinzugefügt), `changed` (Geändert), `fixed` (Behoben), `removed` (Entfernt).
- Inhalt: ein oder mehrere Listenpunkte wie in `CHANGELOG.md` (`- …`, Folgezeilen zwei Leerzeichen
  eingerückt), mit der Issue-Nummer am Ende, z. B. `(#102)`.

Beispiel `changelog.d/102.added.md`:

```markdown
- `rtt-analysis`: OPD bei Austrittspupille im Unendlichen mit Referenzebene (#102).
```

Prüfen und Vorschau:

```bash
python tools/changelog.py check
python tools/changelog.py preview
```

Beim Release führt der Koordinator die Fragmente zusammen und löscht sie:

```bash
python tools/changelog.py release 0.5.0 "GUI-Grundlagen"
```
