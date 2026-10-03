# Agenten einsetzen

Anleitung für die Arbeit mit Programmieragenten (Claude Code) an Raytatouille.
Regeln für die Agenten selbst stehen in [AGENTS.md](../AGENTS.md).

## Rollen

| Rolle | Wer | Aufgabe |
| --- | --- | --- |
| Maintainer | Mensch | entscheidet ADRs, beantwortet Fragen (Label `frage`), merged PRs |
| Koordinator | eine Claude-Code-Sitzung im Haupt-Checkout | Status, freie Issues, PR-Prüfung gegen Abnahmekriterien, neue Issues schneiden; schreibt keinen Produktionscode |
| Implementierer | eine Claude-Code-Sitzung pro Issue, je in eigenem Worktree | Referenztests, Code, PR bis CI grün |
| Reviewer | Subagent `physik-reviewer` (`.claude/agents/`) | prüft jeden Diff vor dem PR |

## Einmalige Einrichtung

1. Claude Code installieren und anmelden: <https://code.claude.com/docs>.
2. GitHub CLI `gh` installieren und anmelden (`gh auth login` oder `GH_TOKEN` setzen).
   Das Token braucht für dieses Repository: Contents, Pull requests, Issues, Workflows
   (Read and write), Actions (Read).
3. Repository klonen und einmal bauen (siehe AGENTS.md, Abschnitt „Bauen und testen“).
4. Im Repository einmal `claude` starten und dem Arbeitsverzeichnis vertrauen.

## Implementierer starten

Pro freiem Issue ein eigenes Terminal:

```bash
claude --worktree issue-<N>
```

Claude Code legt dafür `.claude/worktrees/issue-<N>/` mit eigenem Branch an. Startprompt:

```text
Bearbeite Issue #<N> in PhotonPower/Raytatouille.
1. Lies AGENTS.md, docs/architecture.md und das Issue (gh issue view <N>).
2. Baue im Worktree mit dem Preset dev und prüfe, dass alles grün ist.
3. Schreibe zuerst die Referenztests aus dem Issue, prüfe, dass sie fehlschlagen,
   dann implementiere.
4. Führe alle lokalen Prüfungen aus AGENTS.md aus.
5. Lass den Subagenten physik-reviewer deinen Diff prüfen und behebe alle Pflichtbefunde.
6. Öffne einen PR mit "Closes #<N>" (gh pr create), beobachte die CI
   (gh pr checks --watch) und behebe Fehler, bis alles grün ist.
7. Merge nicht selbst. Melde mir den PR-Link und die Antwort des Reviewers.
Bei offenen Fragen zu Konventionen oder Architektur: Kommentar im Issue, Label frage,
dann anhalten und mich informieren.
```

## Koordinator starten

Im Haupt-Checkout (ohne Worktree):

```bash
claude
```

Startprompt:

```text
Du bist Koordinator für Raytatouille. Du schreibst keinen Produktionscode.
Lies AGENTS.md, docs/agenten.md und docs/architecture.md (Roadmap).
Gib mir für den aktuellen Meilenstein (gh issue list --milestone ..., gh pr list):
- welche Issues frei sind, d. h. alle "Blockiert durch"-Issues sind geschlossen,
- welche PRs offen sind und ob ihre CI grün ist,
- offene Fragen mit Label frage.
Für jeden PR, den ich nenne: Lass physik-reviewer prüfen, gleiche mit den
Abnahmekriterien des Issues ab und gib mir eine Merge-Empfehlung.
Wenn ein Meilenstein abgeschlossen ist, schlage die Issues für den nächsten vor
(Vorlage .github/ISSUE_TEMPLATE/aufgabe.md) und lege sie erst nach meiner Freigabe an.
```

## Ablauf

1. Koordinator nennt freie Issues.
2. Je freies Issue ein Implementierer (`claude --worktree issue-<N>`), höchstens so viele
   wie freie Issues, zu Beginn zwei.
3. Implementierer liefert PR mit grüner CI und Reviewer-Urteil.
4. Koordinator prüft gegen die Abnahmekriterien; du mergst (Squash).
5. `main` ist geschützt: Ein PR muss auf dem aktuellen `main` aufbauen. Ist er veraltet,
   den Implementierer bitten: `git fetch origin && git rebase origin/main`, Tests, Push.
6. Nach dem Merge den Worktree beenden; Claude Code bietet beim Beenden das Aufräumen an.

## Wenn etwas schiefgeht

- **Agent ändert Tests, damit sie grün werden:** PR ablehnen, Regel 2 aus AGENTS.md zitieren.
- **Zwei PRs ändern dieselbe Datei:** den zweiten nach dem Merge des ersten rebasen lassen.
- **Agent rät bei einer Konvention:** Frage-Issue anlegen lassen, Antwort in
  docs/architecture.md oder als ADR festhalten, damit sie für alle gilt.
