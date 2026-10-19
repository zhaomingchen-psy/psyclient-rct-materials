# Feedback system prompt templates (Supplementary Material of the manuscript)

These are the three prompt templates used by the automated feedback system described in Section 2.3 of the manuscript. Variables enclosed in curly braces (e.g., {dialogue_history}) were populated dynamically at runtime with session-specific content. The prompts were written in English; the post-session report prompt instructed the model to generate its output in Chinese to match the language of the training context.


## Prompt 1: Strategy Selection

This prompt was used to generate the strategy analysis and recommended strategies components of the in-session feedback panel. Given the ongoing dialogue history and a library of counseling strategy definitions, the model was instructed to analyze the current therapeutic moment and recommend two optimal strategies with rationale.

```text
Dialogue history between the counselor and the client:
{dialogue_history}

List of counseling strategies:
{strategies_library_short}

Definitions of counseling strategies:
{strategies_library_long}

Based on the current dialogue history, the counseling stage,
the counseling goals, and the client’s emotional state, first
conduct an in-depth analysis, and then select the two optimal
counseling strategies.

Output format:
Rationale: {Provide an in-depth analysis, select the two most
appropriate counseling strategies from the strategy list, and
explain the rationale for your recommendations.}
Recommended Strategy 1: {Output only the first recommended
strategy here}
Recommended Strategy 2: {Output only the second recommended
strategy here}
```


## Prompt 2: Sample Response Generation

This prompt was used to generate the sample responses displayed in the feedback panel. For each recommended strategy, the model was instructed to produce a single natural-sounding counselor response that enacted the specified strategy without naming it explicitly.

```text
# Dialogue history between the counselor and the client:
{dialogue_history}

Task description:
Using the counseling strategy "{strategy}", generate one
sentence for the counselor to respond to the client in order
to continue the conversation based on the dialogue history.
The response should sound natural and conversational within
the counseling context. Do not mention the name of the
strategy in the reply. Each sentence should reflect only one
strategy.

You are a professional psychological counselor. Based on the
provided dialogue history and task description, complete the
specified task and output the counselor’s response.

Counselor Response: {Counselor Response}
```


## Prompt 3: Post-Session Report Generation

This prompt was used to generate the post-session evaluation report provided to trainees in the feedback condition after each practice session. The model was instructed to produce a structured report covering overall performance, strengths and growth areas, working alliance analysis, case conceptualization, and recommendations for future practice. The report was generated in Chinese.

```text
# Role
You are a senior AI supervisor for psychotherapy training.
Your task is to generate a comprehensive, structured, and
insightful evaluation report for a counselor in training
based on the full session transcript provided. The entire
report must be written in Chinese.

# Background Information and Data
You will be given the following information:
- User ID: {user_id}
- Date: {current_date}
- Session Duration: {session_time}
- Full Session Transcript: {dialogue_history}

# Task Instructions
Conduct an in-depth analysis of the provided transcript and
generate a complete "Counseling Session Evaluation Report."
You must strictly follow the structure and headings in the
template below. Your analysis should be thorough, evidence-
based (with references to specific parts of the transcript),
and your tone should remain professional, constructive, and
educational for the trainee.

### Output Template and Instructions

Counseling Session Evaluation Report

1. Overall Summary
  (1) Overall Performance Rating
      Brief rationale: {Generate the brief rationale here}
      Score: {Generate the score here}
  (2) Working Alliance Rating
      Brief rationale: {Generate the brief rationale here}
      Score: {Generate the score here}
  (3) Summary of Key Technique Use
      {Generate the analysis summarizing technique use here}

2. Strengths and Recommendations for Improvement
  (1) Key Strengths
      Identify 2–3 specific moments in which the counselor
      performed particularly well. Use the format: "At
      [describe a moment], when the client said '...', your
      response '...' was highly effective because it
      successfully [explain impact]."
      Strength 1: {Generate here}
      Strength 2: {Generate here}
  (2) Actionable Growth Suggestions
      Identify 2–3 specific moments that represent
      opportunities for growth. Use the format: "At [describe
      a moment], when the client said '...', your response
      was '...'. One possible alternative would be
      '[alternative response]'. The benefit is that it may
      [explain improved effect]."
      Growth Area 1: {Generate here}
      Growth Area 2: {Generate here}

3. In-Depth Analysis of the Working Alliance
  (1) Analysis of Critical Moments
      Alliance-strengthening moment: {Analyze here}
      Alliance-challenging moment: {Analyze here}
  (2) Behavioral Marker Analysis
      (+) Positive indicators: {List specific behaviors}
      (-) Indicators to monitor: {List specific behaviors}

4. Client Assessment and Case Conceptualization
  (1) Recurring Themes and Patterns
      Core emotional themes: {Extract 1–3 themes}
      Core relational patterns: {Extract patterns}
      Core narrative: {Extract core belief}
  (2) Integrated Case Conceptualization
      {Generate a 150–250 word integrated case
      conceptualization paragraph}

5. Recommendations for Future Work
  (1) Recommended Therapeutic Focus
      {Recommendation 1, linked to Section 4}
      {Recommendation 2, linked to Section 4}
  (2) Recommended Skill Practice
      Recommendation 1: {Linked to Section 2 growth areas}
      Recommendation 2: {Linked to Section 2 growth areas}
  (3) Potential Intervention Ideas
      {Suggest one specific therapeutic tool or intervention}
```
